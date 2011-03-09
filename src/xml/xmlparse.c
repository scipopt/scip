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

/**@file   xmldef.h
 * @brief  declarations for XML parsing
 * @author Thorsten Koch
 * @author Marc Pfetsch
 *
 * If SPEC_LIKE_SPACE_HANDLING is not defined, all LF,CR will be changed into spaces and from a
 * sequence of spaces only one will be used.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <blockmemshell/memory.h>

#include "xml.h"
#include "xmldef.h"


#include <sys/types.h>
#ifdef WITH_ZLIB
#include <unistd.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>


#define NAME_EXT_SIZE 128
#define ATTR_EXT_SIZE 4096
#define DATA_EXT_SIZE 4096
#define LINE_BUF_SIZE 8192

#define xml_error(a, b) xml_errmsg(a, b, FALSE, __FILE__, __LINE__)


/* forward declarations */
typedef struct parse_stack_struct PSTACK;
typedef struct parse_pos_struct   PPOS;

/** state of the parser */
enum parse_state_enum
{
   STATE_ERROR,
   STATE_BEFORE,
   STATE_IN_TAG,
   STATE_PCDATA,
   STATE_EOF
};
typedef enum   parse_state_enum   PSTATE;

/** Stack as a (singly) linked list. The top element is the current node. */
struct parse_stack_struct
{
   XML_NODE*             node;
   PSTACK*               next;
};

/** Store the current position in the file and the state of the parser. */
struct parse_pos_struct
{
   const char*           filename;
   FPTYPE                fp;
   char                  buf[LINE_BUF_SIZE];
   int                   pos;
   int                   lineno;
   int                   nextsym;
   int                   lastsym;
   PSTATE                state;
   PSTACK*               top;
};


/** output error message with corresponding line and position */
static void xml_errmsg(
   PPOS*                 ppos,
   const char*           msg,
   XML_Bool              msg_only,
   const char*           file,
   int                   line
   )
{
   assert(ppos != NULL);

   if ( ! msg_only )
   {
      fprintf(stderr, "%s(%d) Error in file %s line %d\n", file, line, ppos->filename, ppos->lineno);

      fprintf(stderr, "%s", ppos->buf);

      if (strchr(ppos->buf, '\n') == NULL)
         fputc('\n', stderr);

      fprintf(stderr, "%*s\n", ppos->pos, "^");
   }
   fprintf(stderr, "%s\n\n", msg);
}


/** Push new element on the parse stack. Ein neues Element auf den Parse Stack legen.
 * TRUE if it worked, FAILURE otherwise.
 */
static
XML_Bool push_pstack(
   PPOS*                 ppos,
   XML_NODE*             node
   )
{
   PSTACK* p;

   assert(ppos != NULL);
   assert(node != NULL);

   debugMessage("Pushing %s\n", node->name);

   ALLOC_FALSE( BMSallocMemory(&p) );
   assert(p != NULL);

   p->node   = node;
   p->next   = ppos->top;
   ppos->top = p;

   return TRUE;
}

/** returns top element on stack (which has to be present) */
static XML_NODE* top_pstack(
   const PPOS*           ppos
   )
{
   assert(ppos      != NULL);
   assert(ppos->top != NULL);

   return ppos->top->node;
}

/** remove top element from stack and deletes it
 *
 * TRUE if ok, FALSE otherwise
 */
static
XML_Bool pop_pstack(
   PPOS*                 ppos                /**< input stream position */
   )
{
   PSTACK* p;
   int     result;

   assert(ppos != NULL);

   if (ppos->top == NULL)
   {
      xml_error(ppos, "Stack underflow");
      result = FALSE;
   }
   else
   {
      result    = TRUE;
      p         = ppos->top;
      ppos->top = p->next;

      debugMessage("Poping %s\n", p->node->name);
      BMSfreeMemory(&p);
   }
   return result;
}

/** remove complete stack */
static
void clear_pstack(
   PPOS*                 ppos
   )
{
   assert(ppos != NULL);

   while (ppos->top != NULL)
      (void)pop_pstack(ppos);
}

/** Returns the next character from the input buffer and fills the buffer if it is empty (similar to
 * fgetc()).
 */
static
int my_getc(
   PPOS*                 ppos
   )
{
   assert(ppos      != NULL);
   assert(ppos->fp  != NULL);
   assert(ppos->pos <  LINE_BUF_SIZE);

   if (ppos->buf[ppos->pos] == '\0')
   {
#if 0
      if (NULL == FGETS(ppos->buf, sizeof(ppos->buf), ppos->fp))
         return EOF;
#else
      int len = FREAD(ppos->buf, sizeof(ppos->buf) - 1, ppos->fp);

      if (len <= 0)
         return EOF;

      ppos->buf[len] = '\0';
#endif
      ppos->pos = 0;
   }
   return (unsigned char)ppos->buf[ppos->pos++];
}


#ifdef SPEC_LIKE_SPACE_HANDLING
/** Read input from fp_in.
 *
 * If there is a LF, CR, CR/LF, or LF/CR it returns exactly on LF.  Also counts the number of
 * characters.
 */
static
int getsymbol(
   PPOS*                 ppos
   )
{
   int c;

   assert(ppos != NULL);

   if (ppos->nextsym == 0)
      c = my_getc(ppos);
   else
   {
      c = ppos->nextsym;
      ppos->nextsym = 0;
   }
   assert(ppos->nextsym == 0);

   if (((c == '\n') && (ppos->lastsym == '\r'))
      || ((c == '\r') && (ppos->lastsym == '\n')))
      c = my_getc(ppos);

   ppos->lastsym = c;

   if (c == '\r')
      c = '\n';

   if (c == '\n')
      ppos->lineno++;

   return c;
}
#else
/** Read input from fp_in (variant).
 *
 *  Here we convert all LF or CR into SPACE and return maximally one SPACE after the other.
 *
 *  @note This function counts lines differently. On systems that have only one '\\r' as line feed
 *  (MAC) it does not count correctly.
 */
static
int getsymbol(
   PPOS*                 ppos
   )
{
   int c;

   assert(ppos != NULL);

   do
   {
      if (ppos->nextsym == 0)
         c = my_getc(ppos);
      else
      {
         c = ppos->nextsym;
         ppos->nextsym = 0;
      }
      assert(ppos->nextsym == 0);

      if (c == '\n')
         ppos->lineno++;

      if ((c == '\n') || (c == '\r'))
         c = ' ';
   } while((c == ' ') && (ppos->lastsym == c));

   ppos->lastsym = c;

   debugMessage("[%c]\n", c);

   return c;
}
#endif

/** Reinserts a character into the input stream */
static
void ungetsymbol(
   PPOS*                 ppos,
   int                   c
   )
{
   assert(ppos          != NULL);
   assert(ppos->nextsym == 0);

   ppos->nextsym = c;
}

/** Skip all spaces and return the next non-space character or EOF */
static
int skip_space(
   PPOS*                 ppos
   )
{
   int c;

   assert(ppos != NULL);

   do { c = getsymbol(ppos); } while(isspace(c));

   return c;
}

/** Get name of a TAG or attribute from the input stream.
 *
 *  Either it returns a pointer to allocated memory which contains the name or it returns NULL if
 *  there is some error.
 */
static
char* get_name(
   PPOS*                 ppos
   )
{
   char*  name = NULL;
   size_t size = 0;
   size_t len  = 0;
   int    c;

   assert(ppos != NULL);

   c = getsymbol(ppos);

   if (!isalpha(c) && (c != '_') && (c != ':'))
   {
      xml_error(ppos, "Name starting with illegal charater");
      return NULL;
   }

   /* The following is wrong: Here almost all characters that we casted to unicode are feasible */
   while (isalnum(c) || (c == '_') || (c == ':') || (c == '.') || (c == '-'))
   {
      if (len + 1 >= size)
      {
         size += NAME_EXT_SIZE;

         if ( name == NULL )
         {
            ALLOC_ABORT( BMSallocMemoryArray(&name, size) );
         }
         else
         {
            ALLOC_ABORT( BMSreallocMemoryArray(&name, size) );
         }
      }
      assert(name != NULL);
      assert(size > len);

      name[len++] = (char)c;

      c = getsymbol(ppos);
   }
   if (c != EOF)
      ungetsymbol(ppos, c);

   assert(name != NULL);

   if (len == 0)
   {
      BMSfreeMemoryArray(&name);
      name = NULL;
   }
   else
      name[len] = '\0';

   return name;
}

/** Read the value of an attribute from the input stream.
 *
 *  The value has to be between two " or ' (the other character is then valid as well). The function
 *  returns a pointer to allocated memory containing the value or it returns NULL in case of an
 *  error.
 */
static
char* get_attrval(
   PPOS*                 ppos
   )
{
   char*  attr = NULL;
   int    c;
   int    stop;
   size_t len = 0;
   size_t size = 0;

   assert(ppos != NULL);

   /* The following is not allowed according to the specification (the value has to be directly
    * after the equation sign). */
   c = skip_space(ppos);

   if ((c != '"') && (c != '\''))
   {
      xml_error(ppos, "Atribute value does not start with \" or \'");
      return NULL;
   }
   stop = c;

   for(;;)
   {
      if (len == size)
      {
         size += ATTR_EXT_SIZE;

         if ( attr == NULL )
         {
            ALLOC_ABORT( BMSallocMemoryArray(&attr, size) );
         }
         else
         {
            ALLOC_ABORT( BMSreallocMemoryArray(&attr, size) );
         }
      }
      assert(attr != NULL);
      assert(size >  len);

      c = getsymbol(ppos);

      if ((c == stop) || (c == EOF))
         break;

      attr[len++] = (char)c;
   }
   if (c != EOF)
      attr[len] = '\0';
   else
   {
      BMSfreeMemoryArray(&attr);
      attr = NULL;
   }
   return attr;
}

/** Skip comment
 *
 *  Return FALSE if an error occurs.
 */
static
int do_comment(
   PPOS*                 ppos
   )
{
   int result = TRUE;
   int c;
   int state = 0;

   assert(ppos != NULL);

   for(;;)
   {
      c = getsymbol(ppos);

      if (c == EOF)
         break;
      if ((c == '>') && (state >= 2))
         break;
      state = (c == '-') ? state + 1 : 0;
   }
   if (c == EOF)
   {
      xml_error(ppos, "Unexpected EOF in comment");
      result = FALSE;
   }
   return result;
}

/** Handles a CDATA section.
 *
 *  Returns a pointer to allocated memory containing the data of this section or NULL in case of an
 *  error.
 */
static
char* do_cdata(
   PPOS*                 ppos
   )
{
   char*  data  = NULL;
   size_t size  = 0;
   size_t len   = 0;
   int    state = 0;
   int    c;

   assert(ppos != NULL);

   for(;;)
   {
      c = getsymbol(ppos);

      if (c == EOF)
         break;

      if (c == ']')
         state++;
      else
         if ((c == '>') && (state >= 2))
            break;
         else
            state = 0;

      if (len == size)
      {
         size += DATA_EXT_SIZE;

         if ( data == NULL )
         {
            ALLOC_ABORT( BMSallocMemoryArray(&data, size) );
         }
         else
         {
            ALLOC_ABORT( BMSreallocMemoryArray(&data, size) );
         }
      }
      assert(data != NULL);
      assert(size >  len);

      data[len++] = (char)c;
   }
   if (c != EOF)
   {
      assert(data != NULL);
      assert(len  >= 2);

      data[len - 2] = '\0';
   }
   else
   {
      assert(data != NULL);

      BMSfreeMemoryArray(&data);
      data = NULL;
      xml_error(ppos, "Unexpected EOF in CDATA");
   }
   return data;
}

/** Handle processing instructions (skipping) */
static
void handle_pi(
   PPOS*                 ppos
   )
{
   int c;

   assert(ppos        != NULL);
   assert(ppos->state == STATE_BEFORE);

   do { c = getsymbol(ppos); } while((c != EOF) && (c != '>'));

   if (c != EOF)
      ppos->state = STATE_PCDATA;
   else
   {
      xml_error(ppos, "Unexpected EOF in PI");
      ppos->state = STATE_ERROR;
   }
}

/** Handles declarations that start with a <!.
 *
 *  This includes comments. Does currenlty not work very well, because of DTDs.
 */
static
void handle_decl(
   PPOS*                 ppos
   )
{
   enum XmlSection {
      IS_COMMENT,
      IS_ATTLIST,
      IS_DOCTYPE,
      IS_ELEMENT,
      IS_ENTITY,
      IS_NOTATION,
      IS_CDATA
   };
   typedef enum XmlSection XMLSECTION;
   
   static struct
   {
      const char* name;
      XMLSECTION  what;
   } key[] =
   {
      { "--",       IS_COMMENT  },
      { "ATTLIST",  IS_ATTLIST  },
      { "DOCTYPE",  IS_DOCTYPE  },
      { "ELEMENT",  IS_ELEMENT  },
      { "ENTITY",   IS_ENTITY   },
      { "NOTATION", IS_NOTATION },
      { "[CDATA[",  IS_CDATA    }
   };
   XML_NODE* node;
   char*   data;
   int     c;
   int     k      = 0;
   int     beg    = 0;
   int     end    = (sizeof(key) / sizeof(key[0])) - 1;

   assert(ppos        != NULL);
   assert(ppos->state == STATE_BEFORE);

   do
   {
      c = getsymbol(ppos);

      for(; (beg <= end) && (c != key[beg].name[k]); beg++)
         ;
      for(; (end >= beg) && (c != key[end].name[k]); end--)
         ;
      k++;
   } while(beg < end);

   if (beg != end)
   {
      xml_error(ppos, "Unknown declaration");

      while((c != EOF) && (c != '>'))
         c = getsymbol(ppos);
   }
   else
   {
      assert(beg == end);
      assert(beg <  (int)(sizeof(key) / sizeof(*key)));

      switch(key[beg].what)
      {
      case IS_COMMENT :
         if (do_comment(ppos))
            ppos->state = STATE_ERROR;
         break;
      case IS_CDATA :
         if ((data = do_cdata(ppos)) == NULL)
            ppos->state = STATE_ERROR;
         else
         {
            if (NULL == (node = xml_new_node("#CDATA", ppos->lineno)))
            {
               xml_error(ppos, "Can't create new node");
               ppos->state = STATE_ERROR;
            }
            else
            {
               BMSduplicateMemoryArray(&node->data, data, strlen(data)+1);
               BMSfreeMemoryArray(&data);
               xml_append_child(top_pstack(ppos), node);
            }
         }
         break;
      case IS_ATTLIST :
      case IS_ELEMENT :
      case IS_NOTATION :
      case IS_ENTITY :
      case IS_DOCTYPE :
         break;
      default :
         abort();
      }
   }
}

/** Handle end tag */
static
void handle_endtag(
   PPOS*                 ppos
   )
{
   char* name;
   int   c;

   assert(ppos != NULL);

   if ((name = get_name(ppos)) == NULL)
      xml_error(ppos, "Missing name in endtag");
   else
   {
      c = skip_space(ppos);

      if (c != '>')
      {
         xml_error(ppos, "Missing '>' in endtag");
         ppos->state = STATE_ERROR;
      }
      else
      {
         if (strcmp(name, top_pstack(ppos)->name))
         {
            xml_error(ppos, "Name of endtag does not match starttag");
            ppos->state = STATE_ERROR;
         }
         else
         {
            if ( pop_pstack(ppos) )
               ppos->state = STATE_PCDATA;
            else
               ppos->state = STATE_ERROR;
         }
         BMSfreeMemoryArray(&name);
      }
   }
}

/** Handle start tag */
static
void handle_starttag(
   PPOS*                 ppos
   )
{
   XML_NODE* node;
   char*   name;

   assert(ppos != NULL);

   if ((name = get_name(ppos)) == NULL)
   {
      xml_error(ppos, "Missing name in tagstart");
      ppos->state = STATE_ERROR;
   }
   else
   {
      if (NULL == (node = xml_new_node(name, ppos->lineno)))
      {
         xml_error(ppos, "Can't create new node");
         ppos->state = STATE_ERROR;
      }
      else
      {
         xml_append_child(top_pstack(ppos), node);

         if ( push_pstack(ppos, node) )
            ppos->state = STATE_IN_TAG;
         else
            ppos->state = STATE_ERROR;
      }
      BMSfreeMemoryArray(&name);
   }
}

/** Checks for next tag */
static
void proc_before(
   PPOS*                 ppos                /**< input stream position */
   )
{
   int c;

   assert(ppos        != NULL);
   assert(ppos->state == STATE_BEFORE);

   c = skip_space(ppos);

   if (c != '<')
   {
      xml_error(ppos, "Expecting '<'");
      ppos->state = STATE_ERROR;
   }
   else
   {
      c = getsymbol(ppos);

      switch(c)
      {
      case EOF :
         xml_error(ppos, "Unexpected EOF");
         ppos->state = STATE_ERROR;
         break;
      case '!' :
         handle_decl(ppos);
         break;
      case '?' :
         handle_pi(ppos);
         break;
      case '/' :
         handle_endtag(ppos);
         break;
      default :
         ungetsymbol(ppos, c);
         handle_starttag(ppos);
         break;
      }
   }
}

/** Process tag */
static
void proc_in_tag(
   PPOS*                 ppos                /**< input stream position */
   )
{
   XML_ATTR* attr;
   int     c;
   int     empty = FALSE;
   char*   name;
   char*   value;

   assert(ppos        != NULL);
   assert(ppos->state == STATE_IN_TAG);

   c = skip_space(ppos);

   if ((c == '/') || (c == '>') || (c == EOF))
   {
      if (c == '/')
      {
         empty = TRUE;
         c = getsymbol(ppos);
      }
      if (c == EOF)
      {
         xml_error(ppos, "Unexpected EOF while in a tag");
         ppos->state = STATE_ERROR;
      }
      if (c == '>')
      {
         ppos->state = STATE_PCDATA;

         if (empty && ! pop_pstack(ppos))
            ppos->state = STATE_ERROR;
      }
      else
      {
         xml_error(ppos, "Expected tag end marker '>'");
         ppos->state = STATE_ERROR;
      }
   }
   else
   {
      ungetsymbol(ppos, c);

      if ((name = get_name(ppos)) == NULL)
      {
         xml_error(ppos, "No name for attribute");
         ppos->state = STATE_ERROR;
      }
      else
      {
         c = skip_space(ppos);

         if ((c != '=') || ((value = get_attrval(ppos)) == NULL))
         {
            xml_error(ppos, "Missing attribute value");
            ppos->state = STATE_ERROR;
            BMSfreeMemoryArray(&name);
         }
         else
         {
            if (NULL == (attr = xml_new_attr(name, value)))
            {
               xml_error(ppos, "Can't create new attribute");
               ppos->state = STATE_ERROR;
            }
            else
            {
               xml_add_attr(top_pstack(ppos), attr);
            }
            BMSfreeMemoryArray(&name);
            BMSfreeMemoryArray(&value);
         }
      }
   }
}

/* Handles PCDATA */
static
void proc_pcdata(
   PPOS*                 ppos                /**< input stream position */
   )
{
   XML_NODE* node;
   char*   data   = NULL;
   size_t  size   = 0;
   size_t  len    = 0;
   int     c;

   assert(ppos        != NULL);
   assert(ppos->state == STATE_PCDATA);

#ifndef SPEC_LIKE_SPACE_HANDLING
   if ((c = skip_space(ppos)) != EOF)
      ungetsymbol(ppos, c);
#endif
   c = getsymbol(ppos);

   while ((c != EOF) && (c != '<'))
   {
      if (len == size - 1) /* leave space for terminating '\0' */
      {
         size += DATA_EXT_SIZE;

         if ( data == NULL )
         {
            ALLOC_ABORT( BMSallocMemoryArray(&data, size) );
         }
         else
         {
            ALLOC_ABORT( BMSreallocMemoryArray(&data, size) );
         }
      }
      assert(data != NULL);
      assert(size > len + 1);

      data[len++] = (char)c;

      c = getsymbol(ppos);
   }
   if (data == NULL)
   {
      if (c == EOF)
         ppos->state = STATE_EOF;
      else if (c == '<')
      {
         ppos->state = STATE_BEFORE;
         ungetsymbol(ppos, c);
      }
      else
      {
         ppos->state = STATE_ERROR;
      }
   }
   else
   {
      assert(len < size);
      data[len] = '\0';

      if (c == EOF)
         ppos->state = STATE_ERROR;
      else
      {
         ungetsymbol(ppos, c);

         if (NULL == (node = xml_new_node("#PCDATA", ppos->lineno)))
         {
            xml_error(ppos, "Can't create new node");
            ppos->state = STATE_ERROR;
         }
         else
         {
            BMSduplicateMemoryArray(&node->data, data, strlen(data)+1);
            xml_append_child(top_pstack(ppos), node);
            ppos->state = STATE_BEFORE;
         }
         BMSfreeMemoryArray(&data);
      }
   }
}

/** Parse input stream */
static
XML_Bool xml_parse(
   PPOS*                 ppos                /**< input stream position */
   )
{
   XML_Bool ok = TRUE;

   while (ok)
   {
      debugMessage("state=%d\n", ppos->state);

      switch (ppos->state)
      {
      case STATE_BEFORE :
         proc_before(ppos);
         break;
      case STATE_IN_TAG :
         proc_in_tag(ppos);
         break;
      case STATE_PCDATA :
         proc_pcdata(ppos);
         break;
      case STATE_EOF :
         ok = FALSE;
         break;
      case STATE_ERROR :
         ok = FALSE;
         break;
      default :
         xml_error(ppos, "Internal Error, illegal state");
         ok = FALSE;
      }
   }
   return (ppos->state == STATE_EOF);
}

/** Parse file */
XML_NODE* xml_process(
   const char*           filename            /**< XML file name */
   )
{
   PPOS      ppos;
   XML_NODE* node = NULL;
   XML_ATTR* attr;
   int       result = FALSE;
   char*     myfilename;

   ALLOC_FALSE( BMSduplicateMemoryArray(&myfilename, filename, strlen(filename) + 5) );

#ifdef WITH_ZLIB
   if (access(filename, R_OK) != 0)
   {
      strcat(myfilename, ".gz");

      /* If .gz also does not work, revert to the old name
       * to get a better error message.
       */
      if (access(myfilename, R_OK) != 0)
         strcpy(myfilename, filename);
   }
#endif
   ppos.fp = FOPEN(myfilename, "r");
   if ( ppos.fp == NULL )
      perror(myfilename);
   else
   {
      ppos.filename = myfilename;
      ppos.buf[0]   = '\0';
      ppos.pos      = 0;
      ppos.lineno   = 1;
      ppos.nextsym  = 0;
      ppos.lastsym  = 0;
      ppos.state    = STATE_BEFORE;
      ppos.top      = NULL;

      node = xml_new_node("#ROOT", ppos.lineno);
      if ( node == NULL )
      {
         xml_error(&ppos, "Can't create new node");
      }
      else
      {
         attr = xml_new_attr("filename", myfilename);
         if ( attr == NULL )
            xml_error(&ppos, "Can't create new attribute");
         else
         {
            xml_add_attr(node, attr);

            /* push root node on stack and start to process */
            if ( push_pstack(&ppos, node) )
            {
               result = xml_parse(&ppos);

               clear_pstack(&ppos);
            }
         }
      }

      if ( ! result && (node != NULL) )
      {
         xml_errmsg(&ppos, "Parsing error, processing stopped", TRUE, __FILE__, __LINE__);
         xml_free_node(node);
         node = NULL;
      }
      if (FCLOSE(ppos.fp))
         perror(myfilename);
   }
   BMSfreeMemoryArray(&myfilename);

   return node;
}






/*----------------------------------------------------------------------------------------------*/


/** create new node */
XML_NODE* xml_new_node(
   const char*           name,
   int                   lineno
   )
{
   XML_NODE* n = NULL;

   assert(name != NULL);

   if ( BMSallocMemory(&n) != NULL )
   {
      BMSclearMemory(n);
      BMSduplicateMemoryArray(&n->name, name, strlen(name)+1);
      n->lineno = lineno;
   }
   return n;
}

/** create new attribute */
XML_ATTR* xml_new_attr(
   const char*           name,
   const char*           value
   )
{
   XML_ATTR* a = NULL;

   assert(name  != NULL);
   assert(value != NULL);

   if ( BMSallocMemory(&a) != NULL )
   {
      BMSclearMemory(a);
      BMSduplicateMemoryArray(&a->name, name, strlen(name)+1);
      BMSduplicateMemoryArray(&a->value, value, strlen(value)+1);
   }
   return a;
}

/** add attribute */
void xml_add_attr(
   XML_NODE*             n,
   XML_ATTR*             a
   )
{
   assert(n != NULL);
   assert(a != NULL);

   a->next      = n->attr_list;
   n->attr_list = a;
}

/** append child node */
void xml_append_child(
   XML_NODE*             parent,
   XML_NODE*             child
   )
{
   assert(parent != NULL);
   assert(child  != NULL);

   child->parent      = parent;
   child->prev_sibl   = parent->last_child;
   child->next_sibl   = NULL;
   parent->last_child = child;

   if (child->prev_sibl != NULL)
      child->prev_sibl->next_sibl = child;

   if (parent->first_child == NULL)
      parent->first_child = child;
}

/** free attribute */
static
void xml_free_attr(
   XML_ATTR*             a
   )
{
   if (a != NULL)
   {
      xml_free_attr(a->next);

      assert(a->name  != NULL);
      assert(a->value != NULL);

      BMSfreeMemoryArray(&a->name);
      BMSfreeMemoryArray(&a->value);
      BMSfreeMemory(&a);
   }
}

/** free node */
void xml_free_node(
   XML_NODE*             n
   )
{
   if (n != NULL)
   {
      xml_free_node(n->first_child);
      xml_free_node(n->next_sibl);
      xml_free_attr(n->attr_list);

      if (n->data != NULL)
      {
         BMSfreeMemoryArray(&n->data);
      }
      assert(n->name != NULL);

      BMSfreeMemoryArray(&n->name);
      BMSfreeMemory(&n);
   }
}

/** output node */
void xml_show_node(
   const XML_NODE*       root
   )
{
   const XML_NODE* n;
   const XML_ATTR* a;

   assert(root != NULL);

   for (n = root; n != NULL; n = n->next_sibl)
   {
      infoMessage("Name: %s\n", n->name);
      infoMessage("Line: %d\n", n->lineno);
      infoMessage("Data: %s\n", (n->data != NULL) ? n->data : "***");

      for(a = n->attr_list; a != NULL; a = a->next)
         infoMessage("Attr: %s = [%s]\n", a->name, a->value);

      if (n->first_child != NULL)
      {
         infoMessage("->\n");
         xml_show_node(n->first_child);
         infoMessage("<-\n");
      }
   }
}

/** get attribute value */
const char* xml_get_attrval(
   const XML_NODE*       node,
   const char*           name
   )
{
   XML_ATTR* a;

   assert(node != NULL);
   assert(name != NULL);

   for (a = node->attr_list; a != NULL; a = a->next)
   {
      if (!strcmp(name, a->name))
         break;
   }

#if 0
   if (a == NULL)
      infoMessage("Error: Attribut %s in TAG <%s> not found\n",
         name, node->name);
#endif

   return (a == NULL) ? NULL : a->value;
}

/** return first node */
const XML_NODE* xml_first_node(
   const XML_NODE*       node,
   const char*           name
   )
{
   const XML_NODE* n;

   assert(node != NULL);
   assert(name != NULL);

   for (n = node; n != NULL; n = n->next_sibl)
   {
      if (!strcmp(name, n->name))
         break;
   }

   return n;
}

/** return next node */
const XML_NODE* xml_next_node(
   const XML_NODE*       node,
   const char*           name
   )
{
   assert(node != NULL);
   assert(name != NULL);

   return (node->next_sibl == NULL) ? NULL : xml_first_node(node->next_sibl, name);
}

/** find next node */
const XML_NODE* xml_find_node(
   const XML_NODE*       node,
   const char*           name
   )
{
   const XML_NODE* n;

   assert(node != NULL);
   assert(name != NULL);

   if (!strcmp(name, node->name))
      return node;

   if (node->first_child != NULL)
   {
      if (NULL != (n = xml_find_node(node->first_child, name)))
         return n;
   }

   if (node->next_sibl != NULL)
   {
      if (NULL != (n = xml_find_node(node->next_sibl, name)))
         return n;
   }

   return NULL;
}

/** return next sibling */
const XML_NODE* xml_next_sibl(
   const XML_NODE*       node
   )
{
   assert(node != NULL);

   return node->next_sibl;
}

/** return previous sibling */
const XML_NODE* xml_prev_sibl(
   const XML_NODE*       node
   )
{
   assert(node != NULL);

   return node->prev_sibl;
}

/** return first child */
const XML_NODE* xml_first_child(
   const XML_NODE*       node
   )
{
   assert(node != NULL);

   return node->first_child;
}

/** return last child */
const XML_NODE* xml_last_child(
   const XML_NODE*       node
   )
{
   assert(node != NULL);

   return node->last_child;
}

/** return name of node */
const char* xml_get_name(
   const XML_NODE*       node
   )
{
   assert(node != NULL);

   return node->name;
}

/** get line number */
int xml_get_line(
   const XML_NODE*       node
   )
{
   assert(node != NULL);

   return node->lineno;
}

/** get data */
const char* xml_get_data(
   const XML_NODE*       node
   )
{
   assert(node != NULL);

   return node->data;
}

/** find PCDATA */
const char* xml_find_pcdata(
   const XML_NODE*       node,
   const char*           name
   )
{
   const XML_NODE* n;

   assert(node != NULL);
   assert(name != NULL);

   if (NULL == (n = xml_find_node(node, name)))
      return NULL;

   if (!strcmp(n->first_child->name, "#PCDATA"))
      return n->first_child->data;

   return NULL;
}

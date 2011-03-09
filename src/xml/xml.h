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

/**@file   xml.h
 * @brief  declarations for XML parsing
 * @author Thorsten Koch
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_XML_H__
#define __SCIP_XML_H__

#ifdef __cplusplus
extern "C" {
#endif


typedef struct XML_ATTR_struct XML_ATTR;

struct XML_ATTR_struct
{
   char*                 name;
   char*                 value;
   XML_ATTR*             next;
};

typedef struct XML_NODE_struct XML_NODE;

struct XML_NODE_struct
{
   char*                 name;
   int                   lineno;
   XML_ATTR*             attr_list;
   XML_NODE*             parent;
   XML_NODE*             prev_sibl;
   XML_NODE*             next_sibl;
   XML_NODE*             first_child;
   XML_NODE*             last_child;
   char*                 data;           /* does not come together with childs */
};

/** Parse file */
extern
XML_NODE* xml_process(
   const char*           filename            /**< XML file name */
   );

/** create new node */
extern
XML_NODE* xml_new_node(
   const char*           name,
   int                   lineno
   );

/** create new attribute */
extern
XML_ATTR* xml_new_attr(
   const char*           name,
   const char*           value
   );

/** add attribute */
extern
void xml_add_attr(
   XML_NODE*             n,
   XML_ATTR*             a
   );

/** append child node */
extern
void xml_append_child(
   XML_NODE*             parent,
   XML_NODE*             child
   );

/** free node */
extern
void xml_free_node(
   XML_NODE*             n
   );

/** output node */
extern
void xml_show_node(
   const XML_NODE*       root
   );

/** get attribute value */
extern
const char* xml_get_attrval(
   const XML_NODE*       node,
   const char*           name
   );

/** return first node */
extern
const XML_NODE* xml_first_node(
   const XML_NODE*       node,
   const char*           name
   );

/** return next node */
extern
const XML_NODE* xml_next_node(
   const XML_NODE*       node,
   const char*           name
   );

/** find next node */
extern
const XML_NODE* xml_find_node(
   const XML_NODE*       node,
   const char*           name
   );

/** return next sibling */
extern
const XML_NODE* xml_next_sibl(
   const XML_NODE*       node
   );

/** return previous sibling */
extern
const XML_NODE* xml_prev_sibl(
   const XML_NODE*       node
   );

/** return first child */
extern
const XML_NODE* xml_first_child(
   const XML_NODE*       node
   );

/** return last child */
extern
const XML_NODE* xml_last_child(
   const XML_NODE*       node
   );

/** return name of node */
extern
const char* xml_get_name(
   const XML_NODE*       node
   );

/** get line number */
extern
int xml_get_line(
   const XML_NODE*       node
   );

/** get data */
extern
const char* xml_get_data(
   const XML_NODE*       node
   );

/** find PCDATA */
extern
const char* xml_find_pcdata(
   const XML_NODE*       node,
   const char*           name
   );

#ifdef __cplusplus
}
#endif

#endif

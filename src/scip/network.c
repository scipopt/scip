/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   network.c
 * @ingroup OTHER_CFILES
 * @brief   Methods for detecting network (sub)matrices
 * @author Rolf van der Hulst
 */

#include "scip/network.h"
#include "scip/scip.h"
#include <assert.h>
#include "blockmemshell/memory.h"

///Types which define matrix sizes
typedef int spqr_matrix_size;
typedef spqr_matrix_size spqr_row;
typedef spqr_matrix_size spqr_col;

#define SPQR_INVALID INT_MAX
#define SPQR_INVALID_ROW SPQR_INVALID
#define SPQR_INVALID_COL SPQR_INVALID


/* Only check in debug mode if the used indices are valid */
#ifndef NDEBUG

/** Determine if the row index is invalid */
static
SCIP_Bool SPQRrowIsInvalid(
   spqr_row              row                 /**< The row to check */
   )
{
   return row == SPQR_INVALID_ROW;
}

/** Determine if the column index is invalid */
static
SCIP_Bool SPQRcolIsInvalid(
   spqr_col              col                 /**< The column to check */
   )
{
   return col == SPQR_INVALID_COL;
}

/** Determine if the row index is valid */
static
SCIP_Bool SPQRrowIsValid(
   spqr_row              row                 /**< The row to check */
   )
{
   return !SPQRrowIsInvalid(row);
}

/** Determine if the column index is valid */
static
SCIP_Bool SPQRcolIsValid(
   spqr_col              col                 /**< The column to check */
   )
{
   return !SPQRcolIsInvalid(col);
}

#endif

/** Columns 0..x correspond to elements 0..x and rows 0..y correspond to elements -1.. -y-1 */
#define MARKER_ROW_ELEMENT (INT_MIN)
#define MARKER_COLUMN_ELEMENT (INT_MAX)
typedef int spqr_element;

/** Checks if an element is a row */
static
SCIP_Bool SPQRelementIsRow(
   spqr_element          element             /**< The element to check */
   )
{
   return element < 0;
}

/** Checks if an element is a column */
static
SCIP_Bool SPQRelementIsColumn(
   spqr_element          element             /**< The element to check */
)
{
   return !SPQRelementIsRow(element);
}

/** Convert an element to the corresponding row index */
static
spqr_row SPQRelementToRow(
   spqr_element          element             /**< The element to convert */
)
{
   assert(SPQRelementIsRow(element));
   return (spqr_row) ( -element - 1 );
}

/** Convert a row to the corresponding element */
static
spqr_element SPQRrowToElement(
   spqr_row              row                 /**< The row to convert */
   )
{
   assert(SPQRrowIsValid(row));
   return (spqr_element) -row - 1;
}

/** Convert an element to the corresponding column index */
static
spqr_col SPQRelementToColumn(
   spqr_element          element             /**< The element to convert */
   )
{
   assert(SPQRelementIsColumn(element));
   return (spqr_col) element;
}

/** Convert a column to the corresponding element */
static
spqr_element SPQRcolumnToElement(
   spqr_col              column              /**< The column to convert */
   )
{
   assert(SPQRcolIsValid(column));
   return (spqr_element) column;
}

/** spqr_member is an index for the members of the SPQR decomposition. The members are the nodes of the SPQR tree.
 * Each member has an associated subgraph, sometimes called a skeleton.
 * If two members are adjacent in the SPQR tree, the two corresponding subgraphs are connected by a 2-separation in any
 * graph represented by the matrix.
 * For members, we reserve all negative values as invalid. We use these negative values in union-find datastructure,
 * where we store the rank of the representative as a negative number.
 */
typedef int spqr_member;
#define SPQR_INVALID_MEMBER (-1)

/** Check if a member is invalid */
static
SCIP_Bool SPQRmemberIsInvalid(
   spqr_member           member              /**< The member to check */
   )
{
   return member < 0;
}

/** Check if a member is valid */
static SCIP_Bool SPQRmemberIsValid(
   spqr_member           member              /**< The member to check */
   )
{
   return !SPQRmemberIsInvalid(member);
}

/** spqr_node is an index for the nodes stored in the decomposition. The nodes are part of each member's skeleton.
 * Similar to spqr_member, we reserve all negative values as invalid and use these in union-find.
 */
typedef int spqr_node;
#define SPQR_INVALID_NODE (-1)

/** Check if a node is invalid */
static
SCIP_Bool SPQRnodeIsInvalid(
   spqr_node             node                /**< The node to check */
)
{
   return node < 0;
}

/** Check if a node is valid */
static
SCIP_Bool SPQRnodeIsValid(
   spqr_node             node                /**< The node to check */
)
{
   return !SPQRnodeIsInvalid(node);
}

/** spqr_arc is an index for the arcs stored in the decomposition. The arcs are part of each member's skeleton.
 * Similar to spqr_node and spqr_member, we reserve all negative values as invalid and use these in some union-find.
 * However, in contrast to spqr_node and spqr_member, the union-find data structure does not represent the arcs,
 * but rather if all arcs that have the same representative we know they are in the same member.
 */
typedef int spqr_arc;
#define SPQR_INVALID_ARC (-1)

/** Check if an arc is invalid */
static
SCIP_Bool SPQRarcIsInvalid(
   spqr_arc              arc                 /**< The arc to check */
   )
{
   return arc < 0;
}

/** Check if an arc is valid */
static
SCIP_Bool SPQRarcIsValid(
   spqr_arc              arc                 /**< The arc to check */
   )
{
   return !SPQRarcIsInvalid(arc);
}

/** The type of the member */
typedef enum
{
   SPQR_MEMBERTYPE_RIGID      = 0,           /** The member's skeleton is 3-connected and has at least 4 edges */
   SPQR_MEMBERTYPE_PARALLEL   = 1,           /** The member's skeleton consists of 2 nodes with at least 3 edges */
   SPQR_MEMBERTYPE_SERIES     = 2,           /** The member's skeleton is a cycle with at least 3 edges */
   SPQR_MEMBERTYPE_LOOP       = 3,           /** The member's skeleton consists of 2 nodes connected by 1 or 2 edges */
   SPQR_MEMBERTYPE_UNASSIGNED = 4            /** To indicate that the member is not a representative member anymore */
} SPQRMemberType;

/** Represents a single node of cyclic doubly-linked list of arc edges*/
typedef struct
{
   spqr_arc previous;
   spqr_arc next;
} SPQRNetworkDecompositionArcListNode;

/** This structure stores the relevant data for a single node. */
typedef struct
{
   spqr_node representativeNode;             /**< Points to the next node in the union-find data structure
                                              * Stores the rank as a negative number if this node represents itself */
   spqr_arc firstArc;                        /**< Points to the head node of the cyclic doubly linked list containing
                                              * the arcs that are adjacent to this node.*/
   int numArcs;                              /**< The number of arcs adjacent to this node */
} SPQRNetworkDecompositionNode;

/** Structure that stores the relevant data for a single arc. */
typedef struct
{
   spqr_node head;                                      /**< The head node of the arc */
   spqr_node tail;                                      /**< The tail node of the arc */
   spqr_member member;                                  /**< The member that contains the arc */
   spqr_member childMember;                             /**< Stores the child member, if this arc points to one */
   SPQRNetworkDecompositionArcListNode headArcListNode; /**< Linked-list node for iterating over the head's arcs */
   SPQRNetworkDecompositionArcListNode tailArcListNode; /**< Linked-list node for iterating over the tail's arcs */
   SPQRNetworkDecompositionArcListNode arcListNode;     /**< Linked-list node for iterating over the member's arcs */

   spqr_element element;                                /**< The element associated to this arc */

   /**Signed union-find for arc directions. If an arc is reversed, it head becomes its tail and vice-versa.
    * For non-rigid members every arc is it's own representative, and the direction is simply given by the boolean.
    * For rigid members, every arc is represented by another arc in the member,
    * and the direction can be found by multiplying the signs along the union-find path
    * We use this data structure to efficiently reverse all arcs in a skeleton.
    */
   spqr_arc representative;                  /**< The representative of the arc */
   SCIP_Bool reversed;                       /**< Whether the arc's head and tail are reversed, or not */
} SPQRNetworkDecompositionArc;

/** Structure that stores the relevant data for a single member */
typedef struct
{
   spqr_member representativeMember;         /** The representative of this member (union-find) */
   SPQRMemberType type;                      /** The type of this member */

   /**The SPQR tree is stored as an arborescence. Each member stores its parents, and each edge of the member
    * pointing to a child member stores the associated member in childMember
    */
   spqr_member parentMember;                 /**< The parent of this member in the arborescence */
   spqr_arc markerToParent;                  /**< The arc pointing to the parent */
   spqr_arc markerOfParent;                  /**< The arc of the parent pointing to this member */

   spqr_arc firstArc;                        /** First arc of the linked list containing the member's arcs */
   int numArcs;                              /** The number of arcs associated to the member */
} SPQRNetworkDecompositionMember;

/** Stores the SPQR forest data structure and its relevant data */
typedef struct
{
   int numArcs;                              /**< The number of slots used in the arc data array */
   int memArcs;                              /**< The amount of space allocated in the arc data array */
   SPQRNetworkDecompositionArc* arcs;        /**< Array of arcs of the SPQR forest, indexed by spqr_arc */
   spqr_arc firstFreeArc;                    /**< Points to the first unused slot in the arcs array */

   int memMembers;                           /**< The amount of space allocated in the member data array */
   int numMembers;                           /**< The number of slots used in the member data array */
   SPQRNetworkDecompositionMember* members;  /**< Array of members of the SPQR forest. Indexed by spqr_member */

   int memNodes;                             /**< The amount of space allocated in the node data array */
   int numNodes;                             /**< The number of slots used in the node data array */
   SPQRNetworkDecompositionNode* nodes;      /**< Array of nodes of the SPQR forest. Indexed by spqr_node */

   int memRows;                              /**< The (maximal) number of rows of the matrix */
   spqr_arc* rowArcs;                        /**< Maps the rows of the matrix to arcs in the decomposition */

   int memColumns;                           /**< The (maximal) number of columns of the matrix */
   spqr_arc* columnArcs;                     /**< Maps the columns of the matrix to arcs in the decomposition */

   SCIP* env;                                /**< SCIP pointer stored for later use */

   int numConnectedComponents;               /** The number of disjoint SPQR trees in the SPQR forest */
} SCIP_NETMATDECDATA;


#ifndef NDEBUG

static SCIP_Bool nodeIsRepresentative(
   const SCIP_NETMATDECDATA* dec,
   spqr_node node
)
{
   assert(dec);
   assert(node < dec->memNodes);
   assert(SPQRnodeIsValid(node));

   return SPQRnodeIsInvalid(dec->nodes[node].representativeNode);
}

#endif

static spqr_node findNode(
   SCIP_NETMATDECDATA* dec,
   spqr_node node
)
{
   assert(dec);
   assert(SPQRnodeIsValid(node));
   assert(node < dec->memNodes);

   spqr_node current = node;
   spqr_node next;

   //traverse down tree to find the root
   while( SPQRnodeIsValid(next = dec->nodes[current].representativeNode))
   {
      current = next;
      assert(current < dec->memNodes);
   }

   spqr_node root = current;
   current = node;

   //update all pointers along path to point to root, flattening the tree
   while( SPQRnodeIsValid(next = dec->nodes[current].representativeNode))
   {
      dec->nodes[current].representativeNode = root;
      current = next;
      assert(current < dec->memNodes);
   }
   return root;
}

static spqr_node findNodeNoCompression(
   const SCIP_NETMATDECDATA* dec,
   spqr_node node
)
{
   assert(dec);
   assert(SPQRnodeIsValid(node));
   assert(node < dec->memNodes);

   spqr_node current = node;
   spqr_node next;

   //traverse down tree to find the root
   while( SPQRnodeIsValid(next = dec->nodes[current].representativeNode))
   {
      current = next;
      assert(current < dec->memNodes);
   }
   spqr_node root = current;
   return root;
}

static spqr_node findArcTail(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);

   spqr_node representative = findNode(dec, dec->arcs[arc].tail);
   dec->arcs[arc].tail = representative;//update the arc information

   return representative;
}

static spqr_node findArcHead(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);

   spqr_node representative = findNode(dec, dec->arcs[arc].head);
   dec->arcs[arc].head = representative;//update the arc information

   return representative;
}

static spqr_node findArcHeadNoCompression(
   const SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);

   spqr_node representative = findNodeNoCompression(dec, dec->arcs[arc].head);
   return representative;
}

static spqr_node findArcTailNoCompression(
   const SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);

   spqr_node representative = findNodeNoCompression(dec, dec->arcs[arc].tail);
   return representative;
}


static spqr_arc getFirstNodeArc(
   const SCIP_NETMATDECDATA* dec,
   spqr_node node
)
{
   assert(dec);
   assert(SPQRnodeIsValid(node));
   assert(node < dec->memNodes);
   return dec->nodes[node].firstArc;
}

static spqr_arc getNextNodeArcNoCompression(
   const SCIP_NETMATDECDATA* dec,
   spqr_arc arc,
   spqr_node node
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);
   assert(nodeIsRepresentative(dec, node));

   if( findArcHeadNoCompression(dec, arc) == node )
   {
      arc = dec->arcs[arc].headArcListNode.next;
   } else
   {
      assert(findArcTailNoCompression(dec, arc) == node);
      arc = dec->arcs[arc].tailArcListNode.next;
   }
   return arc;
}

static spqr_arc getNextNodeArc(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc,
   spqr_node node
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);
   assert(nodeIsRepresentative(dec, node));

   if( findArcHead(dec, arc) == node )
   {
      arc = dec->arcs[arc].headArcListNode.next;
   } else
   {
      assert(findArcTailNoCompression(dec, arc) == node);
      dec->arcs[arc].tail = node;//This assignment is not necessary but speeds up future queries.
      arc = dec->arcs[arc].tailArcListNode.next;
   }
   return arc;
}

static spqr_arc getPreviousNodeArc(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc,
   spqr_node node
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);
   assert(nodeIsRepresentative(dec, node));

   if( findArcHead(dec, arc) == node )
   {
      arc = dec->arcs[arc].headArcListNode.previous;
   } else
   {
      assert(findArcTailNoCompression(dec, arc) == node);
      dec->arcs[arc].tail = node;//This assignment is not necessary but speeds up future queries.
      arc = dec->arcs[arc].tailArcListNode.previous;
   }
   return arc;
}

static void mergeNodeArcList(
   SCIP_NETMATDECDATA* dec,
   spqr_node toMergeInto,
   spqr_node toRemove
)
{

   spqr_arc firstIntoArc = getFirstNodeArc(dec, toMergeInto);
   spqr_arc firstFromArc = getFirstNodeArc(dec, toRemove);
   if( SPQRarcIsInvalid(firstIntoArc))
   {
      //new node has no arcs
      dec->nodes[toMergeInto].numArcs += dec->nodes[toRemove].numArcs;
      dec->nodes[toRemove].numArcs = 0;

      dec->nodes[toMergeInto].firstArc = dec->nodes[toRemove].firstArc;
      dec->nodes[toRemove].firstArc = SPQR_INVALID_ARC;

      return;
   } else if( SPQRarcIsInvalid(firstFromArc))
   {
      //Old node has no arcs; we can just return
      return;
   }

   spqr_arc lastIntoArc = getPreviousNodeArc(dec, firstIntoArc, toMergeInto);
   assert(SPQRarcIsValid(lastIntoArc));
   spqr_arc lastFromArc = getPreviousNodeArc(dec, firstFromArc, toRemove);
   assert(SPQRarcIsValid(lastFromArc));


   SPQRNetworkDecompositionArcListNode* firstIntoNode =
      findArcHead(dec, firstIntoArc) == toMergeInto ? &dec->arcs[firstIntoArc].headArcListNode
                                                    : &dec->arcs[firstIntoArc].tailArcListNode;
   SPQRNetworkDecompositionArcListNode* lastIntoNode =
      findArcHead(dec, lastIntoArc) == toMergeInto ? &dec->arcs[lastIntoArc].headArcListNode
                                                   : &dec->arcs[lastIntoArc].tailArcListNode;

   SPQRNetworkDecompositionArcListNode* firstFromNode =
      findArcHead(dec, firstFromArc) == toRemove ? &dec->arcs[firstFromArc].headArcListNode
                                                 : &dec->arcs[firstFromArc].tailArcListNode;
   SPQRNetworkDecompositionArcListNode* lastFromNode =
      findArcHead(dec, lastFromArc) == toRemove ? &dec->arcs[lastFromArc].headArcListNode
                                                : &dec->arcs[lastFromArc].tailArcListNode;

   firstIntoNode->previous = lastFromArc;
   lastIntoNode->next = firstFromArc;
   firstFromNode->previous = lastIntoArc;
   lastFromNode->next = firstIntoArc;

   dec->nodes[toMergeInto].numArcs += dec->nodes[toRemove].numArcs;
   dec->nodes[toRemove].numArcs = 0;
   dec->nodes[toRemove].firstArc = SPQR_INVALID_ARC;
}

static void arcFlipReversed(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);
   dec->arcs[arc].reversed = !dec->arcs[arc].reversed;
}

static void arcSetReversed(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc,
   SCIP_Bool reversed
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);
   dec->arcs[arc].reversed = reversed;
}

static void arcSetRepresentative(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc,
   spqr_arc representative
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);
   assert(representative == SPQR_INVALID_ARC || SPQRarcIsValid(representative));
   dec->arcs[arc].representative = representative;
}

static spqr_node mergeNodes(
   SCIP_NETMATDECDATA* dec,
   spqr_node first,
   spqr_node second
)
{
   assert(dec);
   assert(nodeIsRepresentative(dec, first));
   assert(nodeIsRepresentative(dec, second));
   assert(first != second);//We cannot merge a node into itself
   assert(first < dec->memNodes);
   assert(second < dec->memNodes);

   //The rank is stored as a negative number: we decrement it making the negative number larger.
   // We want the new root to be the one with 'largest' rank, so smallest number. If they are equal, we decrement.
   spqr_node firstRank = dec->nodes[first].representativeNode;
   spqr_node secondRank = dec->nodes[second].representativeNode;
   if( firstRank > secondRank )
   {
      SCIPswapInts(&first, &second);
   }
   //first becomes representative; we merge all of the arcs of second into first
   mergeNodeArcList(dec, first, second);
   dec->nodes[second].representativeNode = first;
   if( firstRank == secondRank )
   {
      --dec->nodes[first].representativeNode;
   }
   return first;
}

static SCIP_Bool memberIsRepresentative(
   const SCIP_NETMATDECDATA* dec,
   spqr_member member
)
{
   assert(dec);
   assert(member < dec->memMembers);
   assert(SPQRmemberIsValid(member));

   return SPQRmemberIsInvalid(dec->members[member].representativeMember);
}

static spqr_member findMember(
   SCIP_NETMATDECDATA* dec,
   spqr_member member
)
{
   assert(dec);
   assert(member < dec->memMembers);
   assert(SPQRmemberIsValid(member));

   spqr_member current = member;
   spqr_member next;

   //traverse down tree to find the root
   while( SPQRmemberIsValid(next = dec->members[current].representativeMember))
   {
      current = next;
      assert(current < dec->memMembers);
   }

   spqr_member root = current;
   current = member;

   //update all pointers along path to point to root, flattening the tree
   while( SPQRmemberIsValid(next = dec->members[current].representativeMember))
   {
      dec->members[current].representativeMember = root;
      current = next;
      assert(current < dec->memMembers);
   }
   return root;
}

static spqr_member findMemberNoCompression(
   const SCIP_NETMATDECDATA* dec,
   spqr_member member
)
{
   assert(dec);
   assert(member < dec->memMembers);
   assert(SPQRmemberIsValid(member));

   spqr_member current = member;
   spqr_member next;

   //traverse down tree to find the root
   while( SPQRmemberIsValid(next = dec->members[current].representativeMember))
   {
      current = next;
      assert(current < dec->memMembers);
   }

   spqr_member root = current;
   return root;
}

static spqr_member mergeMembers(
   SCIP_NETMATDECDATA* dec,
   spqr_member first,
   spqr_member second
)
{
   assert(dec);
   assert(memberIsRepresentative(dec, first));
   assert(memberIsRepresentative(dec, second));
   assert(first != second);//We cannot merge a member into itself
   assert(first < dec->memMembers);
   assert(second < dec->memMembers);

   //The rank is stored as a negative number: we decrement it making the negative number larger.
   // We want the new root to be the one with 'largest' rank, so smallest number. If they are equal, we decrement.
   spqr_member firstRank = dec->members[first].representativeMember;
   spqr_member secondRank = dec->members[second].representativeMember;
   if( firstRank > secondRank )
   {
      SCIPswapInts(&first, &second);
   }
   dec->members[second].representativeMember = first;
   if( firstRank == secondRank )
   {
      --dec->members[first].representativeMember;
   }
   return first;
}

static spqr_member findArcMember(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);

   spqr_member representative = findMember(dec, dec->arcs[arc].member);
   dec->arcs[arc].member = representative;
   return representative;
}

static spqr_member findArcMemberNoCompression(
   const SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);

   spqr_member representative = findMemberNoCompression(dec, dec->arcs[arc].member);
   return representative;
}

static spqr_member findMemberParent(
   SCIP_NETMATDECDATA* dec,
   spqr_member member
)
{
   assert(dec);
   assert(member < dec->memMembers);
   assert(SPQRmemberIsValid(member));
   assert(memberIsRepresentative(dec, member));


   if( SPQRmemberIsInvalid(dec->members[member].parentMember))
   {
      return dec->members[member].parentMember;
   }
   spqr_member parent_representative = findMember(dec, dec->members[member].parentMember);
   dec->members[member].parentMember = parent_representative;

   return parent_representative;
}

static spqr_member findMemberParentNoCompression(
   const SCIP_NETMATDECDATA* dec,
   spqr_member member
)
{
   assert(dec);
   assert(member < dec->memMembers);
   assert(SPQRmemberIsValid(member));
   assert(memberIsRepresentative(dec, member));

   if( SPQRmemberIsInvalid(dec->members[member].parentMember))
   {
      return dec->members[member].parentMember;
   }
   spqr_member parent_representative = findMemberNoCompression(dec, dec->members[member].parentMember);
   return parent_representative;
}

static spqr_member findArcChildMember(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);

   spqr_member representative = findMember(dec, dec->arcs[arc].childMember);
   dec->arcs[arc].childMember = representative;
   return representative;
}

static spqr_member findArcChildMemberNoCompression(
   const SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);

   spqr_member representative = findMemberNoCompression(dec, dec->arcs[arc].childMember);
   return representative;
}

// Only accounts for CHILD markers, not parent markers!
static SCIP_Bool arcIsMarker(
   const SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);

   return SPQRmemberIsValid(dec->arcs[arc].childMember);
}

static SCIP_Bool arcIsTree(
   const SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);

   return SPQRelementIsRow(dec->arcs[arc].element);
}

typedef struct
{
   spqr_arc representative;
   SCIP_Bool reversed;
} ArcSign;
//find
#ifndef NDEBUG

static SCIP_Bool arcIsRepresentative(
   const SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(arc < dec->memArcs);
   assert(SPQRarcIsValid(arc));

   return SPQRarcIsInvalid(dec->arcs[arc].representative);
}

#endif

static ArcSign findArcSign(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(arc < dec->memArcs);
   assert(SPQRarcIsValid(arc));

   spqr_arc current = arc;
   spqr_arc next;

   SCIP_Bool totalReversed = dec->arcs[current].reversed;
   //traverse down tree to find the root
   while( SPQRarcIsValid(next = dec->arcs[current].representative))
   {
      current = next;
      assert(current < dec->memArcs);
      //swap boolean only if new arc is reversed
      totalReversed = ( totalReversed != dec->arcs[current].reversed );
   }

   spqr_arc root = current;
   current = arc;

   SCIP_Bool currentReversed = totalReversed != dec->arcs[root].reversed;
   //update all pointers along path to point to root, flattening the tree

   while( SPQRarcIsValid(next = dec->arcs[current].representative))
   {
      SCIP_Bool wasReversed = dec->arcs[current].reversed;

      dec->arcs[current].reversed = currentReversed;
      currentReversed = ( currentReversed != wasReversed );

      dec->arcs[current].representative = root;
      current = next;
      assert(current < dec->memArcs);
   }

   ArcSign sign;
   sign.reversed = totalReversed;
   sign.representative = root;
   return sign;
}

static ArcSign findArcSignNoCompression(
   const SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(arc < dec->memArcs);
   assert(SPQRarcIsValid(arc));

   spqr_arc current = arc;
   spqr_arc next;

   SCIP_Bool totalReversed = dec->arcs[current].reversed;
   //traverse down tree to find the root
   while( SPQRarcIsValid(next = dec->arcs[current].representative))
   {
      current = next;
      assert(current < dec->memArcs);
      //swap boolean only if new arc is reversed
      totalReversed = ( totalReversed != dec->arcs[current].reversed );
   }
   ArcSign sign;
   sign.reversed = totalReversed;
   sign.representative = current;
   return sign;
}

//Find the arc tail/head, but accounting for reflection
static spqr_node findEffectiveArcHead(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   if( findArcSign(dec, arc).reversed )
   {
      return findArcTail(dec, arc);
   } else
   {
      return findArcHead(dec, arc);
   }
}

static spqr_node findEffectiveArcTail(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   if( findArcSign(dec, arc).reversed )
   {
      return findArcHead(dec, arc);
   } else
   {
      return findArcTail(dec, arc);
   }
}

static spqr_node findEffectiveArcHeadNoCompression(
   const SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   if( findArcSignNoCompression(dec, arc).reversed )
   {
      return findArcTailNoCompression(dec, arc);
   } else
   {
      return findArcHeadNoCompression(dec, arc);
   }
}

static spqr_node findEffectiveArcTailNoCompression(
   const SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   if( findArcSignNoCompression(dec, arc).reversed )
   {
      return findArcHeadNoCompression(dec, arc);
   } else
   {
      return findArcTailNoCompression(dec, arc);
   }
}

///Merge for signed union find of the arc directions.
///Is not symmetric, in the sense that the arc directions of coponent first are guaranteed not to change but those of second may change
///Based on whether one wants the reflection or not
static spqr_arc mergeArcSigns(
   SCIP_NETMATDECDATA* dec,
   spqr_arc first,
   spqr_arc second,
   SCIP_Bool reflectRelative
)
{
   assert(dec);
   assert(arcIsRepresentative(dec, first));
   assert(arcIsRepresentative(dec, second));
   assert(first != second);//We cannot merge a member into itself
   assert(first < dec->memArcs);
   assert(second < dec->memArcs);

   //The rank is stored as a negative number: we decrement it making the negative number larger.
   // We want the new root to be the one with 'largest' rank, so smallest number. If they are equal, we decrement.
   spqr_member firstRank = dec->arcs[first].representative;
   spqr_member secondRank = dec->arcs[second].representative;

   if( firstRank > secondRank )
   {
      SCIPswapInts(&first, &second);
   }
   dec->arcs[second].representative = first;
   if( firstRank == secondRank )
   {
      --dec->arcs[first].representative;
   }
   //These boolean formula's cover all 16 possible cases, such that the relative orientation of the first is not changed
   SCIP_Bool equal = dec->arcs[first].reversed == dec->arcs[second].reversed;
   dec->arcs[second].reversed = ( equal == reflectRelative );
   if( firstRank > secondRank )
   {
      dec->arcs[first].reversed = ( dec->arcs[first].reversed != reflectRelative );
   }
   return first;
}

static SCIP_Bool arcIsReversedNonRigid(
   const SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);

   return dec->arcs[arc].reversed;
}


static spqr_element arcGetElement(
   const SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);

   return dec->arcs[arc].element;
}

static SCIP_Bool netMatDecDataContainsRow(
   SCIP_NETMATDECDATA*   dec,                /**< The network matrix decomposition */
   int                   row                 /**< The row index that is checked */
)
{
   assert(SPQRrowIsValid(row) && row < dec->memRows);
   assert(dec);
   return SPQRarcIsValid(dec->rowArcs[row]);
}

static SCIP_Bool netMatDecDataContainsColumn(
   SCIP_NETMATDECDATA*   dec,                /**< The network matrix decomposition */
   int                   column              /**< The column index that is checked */
)
{
   assert(SPQRcolIsValid(column) && column < dec->memColumns);
   assert(dec);
   return SPQRarcIsValid(dec->columnArcs[column]);
}

static void setDecompositionColumnArc(
   SCIP_NETMATDECDATA* dec,
   spqr_col col,
   spqr_arc arc
)
{
   assert(SPQRcolIsValid(col) && col < dec->memColumns);
   assert(dec);
   assert(SPQRarcIsValid(arc));
   dec->columnArcs[col] = arc;
}

static void setDecompositionRowArc(
   SCIP_NETMATDECDATA* dec,
   spqr_row row,
   spqr_arc arc
)
{
   assert(SPQRrowIsValid(row) && row < dec->memRows);
   assert(dec);
   assert(SPQRarcIsValid(arc));
   dec->rowArcs[row] = arc;
}

static spqr_arc getDecompositionColumnArc(
   const SCIP_NETMATDECDATA* dec,
   spqr_col col
)
{
   assert(SPQRcolIsValid(col) && col < dec->memColumns);
   assert(dec);
   return dec->columnArcs[col];
}

static spqr_arc getDecompositionRowArc(
   const SCIP_NETMATDECDATA* dec,
   spqr_row row
)
{
   assert(SPQRrowIsValid(row) && row < dec->memRows);
   assert(dec);
   return dec->rowArcs[row];
}

static SCIP_RETCODE netMatDecDataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NETMATDECDATA**  pdec,               /**< buffer to store pointer to created decomposition */
   int                   nrows,              /**< The maximal number of rows that the decomposition can expect */
   int                   ncols               /**< The maximal number of columns that the decomposition can expect */
)
{
   assert(scip);
   assert(pdec);
   assert(!*pdec);

   SCIP_CALL(SCIPallocBlockMemory(scip, pdec));
   SCIP_NETMATDECDATA* dec = *pdec;
   dec->env = scip;

   //Initialize arc array data
   int initialMemArcs = 8;
   {
      assert(initialMemArcs > 0);
      dec->memArcs = initialMemArcs;
      dec->numArcs = 0;
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &dec->arcs, dec->memArcs));
      for( spqr_arc i = 0; i < dec->memArcs; ++i )
      {
         dec->arcs[i].arcListNode.next = i + 1;
         dec->arcs[i].member = SPQR_INVALID_MEMBER;
      }
      dec->arcs[dec->memArcs - 1].arcListNode.next = SPQR_INVALID_ARC;
      dec->firstFreeArc = 0;
   }

   //Initialize member array data
   int initialMemMembers = 8;
   {
      assert(initialMemMembers > 0);
      dec->memMembers = initialMemMembers;
      dec->numMembers = 0;
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &dec->members, dec->memMembers));
   }

   //Initialize node array data
   int initialMemNodes = 8;
   {
      assert(initialMemNodes > 0);
      dec->memNodes = initialMemNodes;
      dec->numNodes = 0;
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &dec->nodes, dec->memNodes));
   }

   //Initialize mappings for rows
   {
      dec->memRows = nrows;
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &dec->rowArcs, dec->memRows));
      for( int i = 0; i < dec->memRows; ++i )
      {
         dec->rowArcs[i] = SPQR_INVALID_ARC;
      }
   }
   //Initialize mappings for columns
   {
      dec->memColumns = ncols;
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &dec->columnArcs, dec->memColumns));
      for( int i = 0; i < dec->memColumns; ++i )
      {
         dec->columnArcs[i] = SPQR_INVALID_ARC;
      }
   }

   dec->numConnectedComponents = 0;
   return SCIP_OKAY;
}

static void netMatDecDataFree(
   SCIP_NETMATDECDATA**      pdec                /**< pointer to the network matrix decomposition to freed */
){
   assert(pdec);
   assert(*pdec);

   SCIP_NETMATDECDATA* dec = *pdec;
   SCIPfreeBlockMemoryArray(dec->env, &dec->columnArcs, dec->memColumns);
   SCIPfreeBlockMemoryArray(dec->env, &dec->rowArcs, dec->memRows);
   SCIPfreeBlockMemoryArray(dec->env, &dec->nodes, dec->memNodes);
   SCIPfreeBlockMemoryArray(dec->env, &dec->members, dec->memMembers);
   SCIPfreeBlockMemoryArray(dec->env, &dec->arcs, dec->memArcs);

   SCIPfreeBlockMemory(dec->env, pdec);
}

static spqr_arc getFirstMemberArc(
   const SCIP_NETMATDECDATA* dec,
   spqr_member member
)
{
   assert(dec);
   assert(SPQRmemberIsValid(member));
   assert(member < dec->memMembers);
   return dec->members[member].firstArc;
}

static spqr_arc getNextMemberArc(
   const SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);
   arc = dec->arcs[arc].arcListNode.next;
   return arc;
}

static spqr_arc getPreviousMemberArc(
   const SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);
   arc = dec->arcs[arc].arcListNode.previous;
   return arc;
}

static void addArcToMemberArcList(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc,
   spqr_member member
)
{
   spqr_arc firstMemberArc = getFirstMemberArc(dec, member);

   if( SPQRarcIsValid(firstMemberArc))
   {
      spqr_arc lastMemberArc = getPreviousMemberArc(dec, firstMemberArc);
      dec->arcs[arc].arcListNode.next = firstMemberArc;
      dec->arcs[arc].arcListNode.previous = lastMemberArc;
      dec->arcs[firstMemberArc].arcListNode.previous = arc;
      dec->arcs[lastMemberArc].arcListNode.next = arc;
   } else
   {
      assert(dec->members[member].numArcs == 0);
      dec->arcs[arc].arcListNode.next = arc;
      dec->arcs[arc].arcListNode.previous = arc;
   }
   dec->members[member].firstArc = arc;
   ++( dec->members[member].numArcs );
}

static SCIP_RETCODE createArc(
   SCIP_NETMATDECDATA* dec,
   spqr_member member,
   SCIP_Bool reversed,
   spqr_arc* pArc
)
{
   assert(dec);
   assert(pArc);
   assert(SPQRmemberIsInvalid(member) || memberIsRepresentative(dec, member));

   spqr_arc idx = dec->firstFreeArc;
   if( SPQRarcIsValid(idx))
   {
      dec->firstFreeArc = dec->arcs[idx].arcListNode.next;
   } else
   {
      //Enlarge array, no free nodes in arc list
      int newSize = 2 * dec->memArcs;
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &dec->arcs, dec->memArcs, newSize));
      for( int i = dec->memArcs + 1; i < newSize; ++i )
      {
         dec->arcs[i].arcListNode.next = i + 1;
         dec->arcs[i].member = SPQR_INVALID_MEMBER;
      }
      dec->arcs[newSize - 1].arcListNode.next = SPQR_INVALID_ARC;
      dec->firstFreeArc = dec->memArcs + 1;
      idx = dec->memArcs;
      dec->memArcs = newSize;
   }

   dec->arcs[idx].tail = SPQR_INVALID_NODE;
   dec->arcs[idx].head = SPQR_INVALID_NODE;
   dec->arcs[idx].member = member;
   dec->arcs[idx].childMember = SPQR_INVALID_MEMBER;
   dec->arcs[idx].reversed = reversed;

   dec->arcs[idx].headArcListNode.next = SPQR_INVALID_ARC;
   dec->arcs[idx].headArcListNode.previous = SPQR_INVALID_ARC;
   dec->arcs[idx].tailArcListNode.next = SPQR_INVALID_ARC;
   dec->arcs[idx].tailArcListNode.previous = SPQR_INVALID_ARC;

   dec->numArcs++;

   *pArc = idx;

   return SCIP_OKAY;
}

static SCIP_RETCODE createRowArc(
   SCIP_NETMATDECDATA* dec,
   spqr_member member,
   spqr_arc* pArc,
   spqr_row row,
   SCIP_Bool reversed
)
{
   SCIP_CALL(createArc(dec, member, reversed, pArc));
   setDecompositionRowArc(dec, row, *pArc);
   addArcToMemberArcList(dec, *pArc, member);
   dec->arcs[*pArc].element = SPQRrowToElement(row);

   return SCIP_OKAY;
}

static SCIP_RETCODE createColumnArc(
   SCIP_NETMATDECDATA* dec,
   spqr_member member,
   spqr_arc* pArc,
   spqr_col column,
   SCIP_Bool reversed
)
{
   SCIP_CALL(createArc(dec, member, reversed, pArc));
   setDecompositionColumnArc(dec, column, *pArc);
   addArcToMemberArcList(dec, *pArc, member);
   dec->arcs[*pArc].element = SPQRcolumnToElement(column);

   return SCIP_OKAY;
}

static SCIP_RETCODE createMember(
   SCIP_NETMATDECDATA* dec,
   SPQRMemberType type,
   spqr_member* pMember
)
{
   assert(dec);
   assert(pMember);

   if( dec->numMembers == dec->memMembers )
   {
      int newSize = dec->memMembers * 2;
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &dec->members, dec->memMembers, newSize));
      dec->memMembers = newSize;
   }
   SPQRNetworkDecompositionMember* data = &dec->members[dec->numMembers];
   data->markerOfParent = SPQR_INVALID_ARC;
   data->markerToParent = SPQR_INVALID_ARC;
   data->firstArc = SPQR_INVALID_ARC;
   data->representativeMember = SPQR_INVALID_MEMBER;
   data->numArcs = 0;
   data->parentMember = SPQR_INVALID_MEMBER;
   data->type = type;

   *pMember = dec->numMembers;

   dec->numMembers++;
   return SCIP_OKAY;
}

static SCIP_RETCODE createNode(
   SCIP_NETMATDECDATA* dec,
   spqr_node* pNode
)
{

   if( dec->numNodes == dec->memNodes )
   {
      int newSize = dec->memNodes * 2;
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &dec->nodes, dec->memNodes, newSize));
      dec->memNodes = newSize;
   }
   *pNode = dec->numNodes;
   dec->nodes[dec->numNodes].representativeNode = SPQR_INVALID_NODE;
   dec->nodes[dec->numNodes].firstArc = SPQR_INVALID_ARC;
   dec->nodes[dec->numNodes].numArcs = 0;
   dec->numNodes++;

   return SCIP_OKAY;
}

static void removeArcFromNodeArcList(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc,
   spqr_node node,
   SCIP_Bool nodeIsHead
)
{
   SPQRNetworkDecompositionArcListNode* arcListNode = nodeIsHead ? &dec->arcs[arc].headArcListNode
                                                                 : &dec->arcs[arc].tailArcListNode;

   if( dec->nodes[node].numArcs == 1 )
   {
      dec->nodes[node].firstArc = SPQR_INVALID_ARC;
   } else
   {
      spqr_arc next_arc = arcListNode->next;
      spqr_arc prev_arc = arcListNode->previous;
      SPQRNetworkDecompositionArcListNode* nextListNode =
         findArcHead(dec, next_arc) == node ? &dec->arcs[next_arc].headArcListNode
                                            : &dec->arcs[next_arc].tailArcListNode;
      SPQRNetworkDecompositionArcListNode* prevListNode =
         findArcHead(dec, prev_arc) == node ? &dec->arcs[prev_arc].headArcListNode
                                            : &dec->arcs[prev_arc].tailArcListNode;

      nextListNode->previous = prev_arc;
      prevListNode->next = next_arc;

      if( dec->nodes[node].firstArc == arc )
      {
         dec->nodes[node].firstArc = next_arc;
      }
   }
   --( dec->nodes[node].numArcs );
}

static void addArcToNodeArcList(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc,
   spqr_node node,
   SCIP_Bool nodeIsHead
)
{
   assert(nodeIsRepresentative(dec, node));

   spqr_arc firstNodeArc = getFirstNodeArc(dec, node);

   SPQRNetworkDecompositionArcListNode* arcListNode = nodeIsHead ? &dec->arcs[arc].headArcListNode
                                                                 : &dec->arcs[arc].tailArcListNode;
   if( SPQRarcIsValid(firstNodeArc))
   {
      SCIP_Bool nextIsHead = findArcHead(dec, firstNodeArc) == node;
      SPQRNetworkDecompositionArcListNode* nextListNode = nextIsHead ? &dec->arcs[firstNodeArc].headArcListNode
                                                                     : &dec->arcs[firstNodeArc].tailArcListNode;
      spqr_arc lastNodeArc = nextListNode->previous;

      arcListNode->next = firstNodeArc;
      arcListNode->previous = lastNodeArc;


      SCIP_Bool previousIsHead = findArcHead(dec, lastNodeArc) == node;
      SPQRNetworkDecompositionArcListNode* previousListNode = previousIsHead ? &dec->arcs[lastNodeArc].headArcListNode
                                                                             : &dec->arcs[lastNodeArc].tailArcListNode;
      previousListNode->next = arc;
      nextListNode->previous = arc;

   } else
   {
      arcListNode->next = arc;
      arcListNode->previous = arc;
   }
   dec->nodes[node].firstArc = arc;
   ++dec->nodes[node].numArcs;
   if( nodeIsHead )
   {
      dec->arcs[arc].head = node;
   } else
   {
      dec->arcs[arc].tail = node;
   }
}

static void setArcHeadAndTail(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc,
   spqr_node head,
   spqr_node tail
)
{
   addArcToNodeArcList(dec, arc, head, TRUE);
   addArcToNodeArcList(dec, arc, tail, FALSE);
}

static void clearArcHeadAndTail(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc
)
{
   removeArcFromNodeArcList(dec, arc, findArcHead(dec, arc), TRUE);
   removeArcFromNodeArcList(dec, arc, findArcTail(dec, arc), FALSE);
   dec->arcs[arc].head = SPQR_INVALID_NODE;
   dec->arcs[arc].tail = SPQR_INVALID_NODE;
}

static void changeArcHead(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc,
   spqr_node oldHead,
   spqr_node newHead
)
{
   assert(nodeIsRepresentative(dec, oldHead));
   assert(nodeIsRepresentative(dec, newHead));
   removeArcFromNodeArcList(dec, arc, oldHead, TRUE);
   addArcToNodeArcList(dec, arc, newHead, TRUE);
}

static void changeArcTail(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc,
   spqr_node oldTail,
   spqr_node newTail
)
{
   assert(nodeIsRepresentative(dec, oldTail));
   assert(nodeIsRepresentative(dec, newTail));
   removeArcFromNodeArcList(dec, arc, oldTail, FALSE);
   addArcToNodeArcList(dec, arc, newTail, FALSE);
}

static int nodeDegree(
   SCIP_NETMATDECDATA* dec,
   spqr_node node
)
{
   assert(dec);
   assert(SPQRnodeIsValid(node));
   assert(node < dec->memNodes);
   return dec->nodes[node].numArcs;
}

static SPQRMemberType getMemberType(
   const SCIP_NETMATDECDATA* dec,
   spqr_member member
)
{
   assert(dec);
   assert(SPQRmemberIsValid(member));
   assert(member < dec->memMembers);
   assert(memberIsRepresentative(dec, member));
   return dec->members[member].type;
}

static void updateMemberType(
   const SCIP_NETMATDECDATA* dec,
   spqr_member member,
   SPQRMemberType type
)
{
   assert(dec);
   assert(SPQRmemberIsValid(member));
   assert(member < dec->memMembers);
   assert(memberIsRepresentative(dec, member));

   dec->members[member].type = type;
}

static spqr_arc markerToParent(
   const SCIP_NETMATDECDATA* dec,
   spqr_member member
)
{
   assert(dec);
   assert(SPQRmemberIsValid(member));
   assert(member < dec->memMembers);
   assert(memberIsRepresentative(dec, member));
   return dec->members[member].markerToParent;
}

static void updateMemberParentInformation(
   SCIP_NETMATDECDATA* dec,
   const spqr_member newMember,
   const spqr_member toRemove
)
{
   assert(memberIsRepresentative(dec, newMember));
   assert(findMemberNoCompression(dec, toRemove) == newMember);

   dec->members[newMember].markerOfParent = dec->members[toRemove].markerOfParent;
   dec->members[newMember].markerToParent = dec->members[toRemove].markerToParent;
   dec->members[newMember].parentMember = dec->members[toRemove].parentMember;

   dec->members[toRemove].markerOfParent = SPQR_INVALID_ARC;
   dec->members[toRemove].markerToParent = SPQR_INVALID_ARC;
   dec->members[toRemove].parentMember = SPQR_INVALID_MEMBER;
}

static spqr_arc markerOfParent(
   const SCIP_NETMATDECDATA* dec,
   spqr_member member
)
{
   assert(dec);
   assert(SPQRmemberIsValid(member));
   assert(member < dec->memMembers);
   assert(memberIsRepresentative(dec, member));
   return dec->members[member].markerOfParent;
}


static int getNumMemberArcs(
   const SCIP_NETMATDECDATA* dec,
   spqr_member member
)
{
   assert(dec);
   assert(SPQRmemberIsValid(member));
   assert(member < dec->memMembers);
   assert(memberIsRepresentative(dec, member));
   return dec->members[member].numArcs;
}

static int getNumNodes(const SCIP_NETMATDECDATA* dec)
{
   assert(dec);
   return dec->numNodes;
}

static int getNumMembers(const SCIP_NETMATDECDATA* dec)
{
   assert(dec);
   return dec->numMembers;
}

static SCIP_RETCODE createStandaloneParallel(
   SCIP_NETMATDECDATA* dec,
   spqr_col* columns,
   SCIP_Bool* reversed,
   int num_columns,
   spqr_row row,
   spqr_member* pMember
)
{
   spqr_member member;
   SCIP_CALL(createMember(dec, num_columns <= 1 ? SPQR_MEMBERTYPE_LOOP : SPQR_MEMBERTYPE_PARALLEL, &member));

   spqr_arc row_arc;
   SCIP_CALL(createRowArc(dec, member, &row_arc, row, num_columns <= 1));

   spqr_arc col_arc;
   for( int i = 0; i < num_columns; ++i )
   {
      SCIP_CALL(createColumnArc(dec, member, &col_arc, columns[i], reversed[i]));
   }
   *pMember = member;

   ++dec->numConnectedComponents;
   return SCIP_OKAY;
}

//TODO: fix tracking connectivity more cleanly, should not be left up to the algorithms ideally
static SCIP_RETCODE createConnectedParallel(
   SCIP_NETMATDECDATA* dec,
   spqr_col* columns,
   SCIP_Bool* reversed,
   int num_columns,
   spqr_row row,
   spqr_member* pMember
)
{
   spqr_member member;
   SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &member));

   spqr_arc row_arc;
   SCIP_CALL(createRowArc(dec, member, &row_arc, row, FALSE));

   spqr_arc col_arc;
   for( int i = 0; i < num_columns; ++i )
   {
      SCIP_CALL(createColumnArc(dec, member, &col_arc, columns[i], reversed[i]));
   }
   *pMember = member;

   return SCIP_OKAY;
}

static SCIP_RETCODE createStandaloneSeries(
   SCIP_NETMATDECDATA* dec,
   spqr_row* rows,
   SCIP_Bool* reversed,
   int numRows,
   spqr_col col,
   spqr_member* pMember
)
{
   spqr_member member;
   SCIP_CALL(createMember(dec, numRows <= 1 ? SPQR_MEMBERTYPE_LOOP : SPQR_MEMBERTYPE_SERIES, &member));

   spqr_arc colArc;
   SCIP_CALL(createColumnArc(dec, member, &colArc, col, FALSE));

   spqr_arc rowArc;
   for( int i = 0; i < numRows; ++i )
   {
      SCIP_CALL(createRowArc(dec, member, &rowArc, rows[i], !reversed[i]));
   }
   *pMember = member;
   ++dec->numConnectedComponents;
   return SCIP_OKAY;
}

static SCIP_RETCODE createConnectedSeries(
   SCIP_NETMATDECDATA* dec,
   spqr_row* rows,
   SCIP_Bool* reversed,
   int numRows,
   spqr_col col,
   spqr_member* pMember
)
{
   spqr_member member;
   SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &member));

   spqr_arc colArc;
   SCIP_CALL(createColumnArc(dec, member, &colArc, col, FALSE));

   spqr_arc rowArc;
   for( int i = 0; i < numRows; ++i )
   {
      SCIP_CALL(createRowArc(dec, member, &rowArc, rows[i], !reversed[i]));
   }
   *pMember = member;
   return SCIP_OKAY;
}

static void removeArcFromMemberArcList(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc,
   spqr_member member
)
{
   assert(findArcMemberNoCompression(dec, arc) == member);
   assert(memberIsRepresentative(dec, member));

   if( dec->members[member].numArcs == 1 )
   {
      dec->members[member].firstArc = SPQR_INVALID_ARC;

   } else
   {
      spqr_arc nextArc = dec->arcs[arc].arcListNode.next;
      spqr_arc prevArc = dec->arcs[arc].arcListNode.previous;

      dec->arcs[nextArc].arcListNode.previous = prevArc;
      dec->arcs[prevArc].arcListNode.next = nextArc;

      if( dec->members[member].firstArc == arc )
      {
         dec->members[member].firstArc = nextArc;
      }
   }


   --( dec->members[member].numArcs );
}

typedef struct
{
   spqr_arc arc;
   SCIP_Bool reversed;
} FindCycleCall;

static void process_arc(
   spqr_row* fundamental_cycle_arcs,
   int* num_cycle_arcs,
   FindCycleCall* callStack,
   int* callStackSize,
   spqr_arc arc,
   const SCIP_NETMATDECDATA* dec,
   SCIP_Bool* fundamental_cycle_direction,
   SCIP_Bool arcIsReversed
)
{
   assert(arcIsTree(dec, arc));
   if( !arcIsMarker(dec, arc))
   {
      spqr_member current_member = findArcMemberNoCompression(dec, arc);
      if( markerToParent(dec, current_member) == arc )
      {
         spqr_arc other_arc = markerOfParent(dec, current_member);
         assert(!arcIsTree(dec, other_arc));
         callStack[*callStackSize].arc = other_arc;
         callStack[*callStackSize].reversed = arcIsReversed;
         ++( *callStackSize );
      } else
      {
         spqr_element element = arcGetElement(dec, arc);
         assert(SPQRelementIsRow(element));
         spqr_row row = SPQRelementToRow(element);
         fundamental_cycle_arcs[*num_cycle_arcs] = row;
         fundamental_cycle_direction[*num_cycle_arcs] = arcIsReversed;
         ++( *num_cycle_arcs );
      }
   } else
   {
      spqr_member child_member = findArcChildMemberNoCompression(dec, arc);
      spqr_arc other_arc = markerToParent(dec, child_member);
      assert(!arcIsTree(dec, other_arc));
      callStack[*callStackSize].arc = other_arc;
      callStack[*callStackSize].reversed = arcIsReversed;
      ++( *callStackSize );
   }
}

static int decompositionGetFundamentalCycleRows(
   const SCIP_NETMATDECDATA* dec,
   spqr_col column,
   spqr_row* output,
   SCIP_Bool* computedSignStorage
)
{
   spqr_arc arc = getDecompositionColumnArc(dec, column);
   if( SPQRarcIsInvalid(arc))
   {
      return 0;
   }
   int num_rows = 0;

   FindCycleCall* callStack;
   SCIP_RETCODE result = SCIPallocBlockMemoryArray(dec->env, &callStack, dec->memRows);
   if( result != SCIP_OKAY )
   {
      return -1;
   }
   int callStackSize = 1;
   callStack[0].arc = arc;
   callStack[0].reversed = FALSE;

   SCIP_Bool* nodeVisited;
   result = SCIPallocBlockMemoryArray(dec->env, &nodeVisited, dec->numNodes);
   if( result != SCIP_OKAY )
   {
      return -1;
   }
   for( int i = 0; i < dec->numNodes; ++i )
   {
      nodeVisited[i] = FALSE;
   }

   typedef struct
   {
      spqr_node node;
      spqr_arc nodeArc;
   } DFSCallData;
   DFSCallData* pathSearchCallStack;
   result = SCIPallocBlockMemoryArray(dec->env, &pathSearchCallStack, dec->numNodes);
   if( result != SCIP_OKAY )
   {
      return -1;
   }
   int pathSearchCallStackSize = 0;

   while( callStackSize > 0 )
   {
      spqr_arc column_arc = callStack[callStackSize - 1].arc;
      SCIP_Bool reverseEverything = callStack[callStackSize - 1].reversed;
      --callStackSize;
      spqr_member column_arc_member = findArcMemberNoCompression(dec, column_arc);
      switch( getMemberType(dec, column_arc_member))
      {
         case SPQR_MEMBERTYPE_RIGID:
         {
            spqr_node source = findEffectiveArcTailNoCompression(dec, column_arc);
            spqr_node target = findEffectiveArcHeadNoCompression(dec, column_arc);

            assert(pathSearchCallStackSize == 0);
            pathSearchCallStack[0].node = source;
            pathSearchCallStack[0].nodeArc = getFirstNodeArc(dec, source);
            pathSearchCallStackSize++;
            while( pathSearchCallStackSize > 0 )
            {
               DFSCallData* dfsData = &pathSearchCallStack[pathSearchCallStackSize - 1];
               nodeVisited[dfsData->node] = TRUE;
               //cannot be a tree arc which is its parent
               if( arcIsTree(dec, dfsData->nodeArc) &&
                   ( pathSearchCallStackSize <= 1 ||
                     dfsData->nodeArc != pathSearchCallStack[pathSearchCallStackSize - 2].nodeArc ))
               {
                  spqr_node head = findEffectiveArcHeadNoCompression(dec, dfsData->nodeArc);
                  spqr_node tail = findEffectiveArcTailNoCompression(dec, dfsData->nodeArc);
                  spqr_node other = head == dfsData->node ? tail : head;
                  assert(other != dfsData->node);
                  assert(!nodeVisited[other]);
                  if( other == target )
                  {
                     break;
                  }
                  //We go up a level: add new node to the call stack

                  pathSearchCallStack[pathSearchCallStackSize].node = other;
                  pathSearchCallStack[pathSearchCallStackSize].nodeArc = getFirstNodeArc(dec, other);
                  ++pathSearchCallStackSize;
                  continue;
               }
               do
               {
                  dfsData->nodeArc = getNextNodeArcNoCompression(dec, dfsData->nodeArc, dfsData->node);
                  if( dfsData->nodeArc == getFirstNodeArc(dec, dfsData->node))
                  {
                     --pathSearchCallStackSize;
                     dfsData = &pathSearchCallStack[pathSearchCallStackSize - 1];
                  } else
                  {
                     break;
                  }
               } while( pathSearchCallStackSize > 0 );
            }
            for( int i = 0; i < pathSearchCallStackSize; ++i )
            {
               if( arcIsTree(dec, pathSearchCallStack[i].nodeArc))
               {
                  SCIP_Bool arcReversedInPath =
                     findEffectiveArcHeadNoCompression(dec, pathSearchCallStack[i].nodeArc) ==
                     pathSearchCallStack[i].node;
                  process_arc(output, &num_rows, callStack, &callStackSize, pathSearchCallStack[i].nodeArc, dec,
                              computedSignStorage, arcReversedInPath != reverseEverything);
               }
            }

            pathSearchCallStackSize = 0;
            break;
         }
         case SPQR_MEMBERTYPE_PARALLEL:
         {
            SCIP_Bool columnReversed = arcIsReversedNonRigid(dec, column_arc);

            spqr_arc first_arc = getFirstMemberArc(dec, column_arc_member);
            spqr_arc iter_arc = first_arc;
            int tree_count = 0;
            do
            {
               if( arcIsTree(dec, iter_arc))
               {
                  SCIP_Bool treeIsReversed = arcIsReversedNonRigid(dec, iter_arc);
                  process_arc(output, &num_rows, callStack, &callStackSize, iter_arc, dec,
                              computedSignStorage, ( columnReversed != treeIsReversed ) != reverseEverything);
                  tree_count++;
               }
               iter_arc = getNextMemberArc(dec, iter_arc);
            } while( iter_arc != first_arc );
            if( tree_count != 1 )
            {
               return -1;
            }
            break;
         }
         case SPQR_MEMBERTYPE_LOOP:
         case SPQR_MEMBERTYPE_SERIES:
         {
            SCIP_Bool columnReversed = arcIsReversedNonRigid(dec, column_arc);
            spqr_arc first_arc = getFirstMemberArc(dec, column_arc_member);
            spqr_arc iter_arc = first_arc;
            int nontree_count = 0;
            do
            {
               if( arcIsTree(dec, iter_arc))
               {
                  SCIP_Bool treeIsReversed = arcIsReversedNonRigid(dec, iter_arc);
                  process_arc(output, &num_rows, callStack, &callStackSize, iter_arc, dec,
                              computedSignStorage, ( columnReversed == treeIsReversed ) != reverseEverything);
               } else
               {
                  nontree_count++;
               }
               iter_arc = getNextMemberArc(dec, iter_arc);
            } while( iter_arc != first_arc );
            if( nontree_count != 1 )
            {
               return -1;
            }
            break;
         }
         case SPQR_MEMBERTYPE_UNASSIGNED:
            assert(FALSE);
      }
   }
   SCIPfreeBlockMemoryArray(dec->env, &pathSearchCallStack, dec->numNodes);
   SCIPfreeBlockMemoryArray(dec->env, &nodeVisited, dec->numNodes);
   SCIPfreeBlockMemoryArray(dec->env, &callStack, dec->memRows);
   return num_rows;
}

static SCIP_Bool netMatDecDataVerifyCycle(
   SCIP* scip,
   const SCIP_NETMATDECDATA* dec,
   int column,
   const int* nonzrowidx,
   const double* nonzvals,
   int num_rows,
   int* pathrowstorage,
   SCIP_Bool* pathsignstorage
)
{
   int num_found_rows = decompositionGetFundamentalCycleRows(dec, column, pathrowstorage, pathsignstorage);

   if( num_found_rows != num_rows )
   {
      return FALSE;
   }
   if( num_rows == 0 )
   {
      return TRUE;
   }
   spqr_row * pathRow;
   int * pathRowReversed;

   SCIP_RETCODE code = SCIPallocBufferArray(scip, &pathRow, num_rows);
   if( code != SCIP_OKAY )
   {
      return FALSE;
   }
   code = SCIPallocBufferArray(scip, &pathRowReversed, num_rows);
   if( code != SCIP_OKAY )
   {
      SCIPfreeBufferArray(scip,&pathRow);
      return FALSE;
   }
   for( int i = 0; i < num_rows; ++i )
   {
      pathRow[i] = pathrowstorage[i];
      pathRowReversed[i] = pathsignstorage[i] ? 1 : 0;
   }
   SCIPsortIntInt(pathRow,pathRowReversed,num_rows);

   spqr_row * secondPathRow;
   int * secondPathRowReversed;
   code = SCIPallocBufferArray(scip, &secondPathRow, num_rows);
   if( code != SCIP_OKAY )
   {
      SCIPfreeBufferArray(scip,&pathRow);
      SCIPfreeBufferArray(scip,&pathRowReversed);
      return FALSE;
   }
   code = SCIPallocBufferArray(scip, &secondPathRowReversed, num_rows);
   if( code != SCIP_OKAY )
   {
      SCIPfreeBufferArray(scip,&pathRow);
      SCIPfreeBufferArray(scip,&pathRowReversed);
      SCIPfreeBufferArray(scip,&secondPathRow);
      return FALSE;
   }
   for( int i = 0; i < num_rows; ++i )
   {
      secondPathRow[i] = nonzrowidx[i];
      secondPathRowReversed[i] = nonzvals[i] < 0.0 ? 1 : 0;
   }

   SCIPsortIntInt(secondPathRow,secondPathRowReversed,num_rows);

   SCIP_Bool good = TRUE;
   for( int i = 0; i < num_rows; ++i )
   {
      if( pathRow[i] != secondPathRow[i] || pathRowReversed[i] != secondPathRowReversed[i] )
      {
         good = FALSE;
         break;
      }
   }
   SCIPfreeBufferArray(scip, &secondPathRowReversed);
   SCIPfreeBufferArray(scip, &secondPathRow);
   SCIPfreeBufferArray(scip, &pathRowReversed);
   SCIPfreeBufferArray(scip, &pathRow);
   return good;
}

static spqr_member largestMemberID(const SCIP_NETMATDECDATA* dec)
{
   return dec->numMembers;
}

static spqr_arc largestArcID(const SCIP_NETMATDECDATA* dec)
{
   return dec->numArcs;
}

static spqr_node largestNodeID(const SCIP_NETMATDECDATA* dec)
{
   return dec->numNodes;
}

static int numConnectedComponents(const SCIP_NETMATDECDATA* dec)
{
   return dec->numConnectedComponents;
}

static SCIP_RETCODE createChildMarker(
   SCIP_NETMATDECDATA* dec,
   spqr_member member,
   spqr_member child,
   SCIP_Bool isTree,
   spqr_arc* pArc,
   SCIP_Bool reversed
)
{
   SCIP_CALL(createArc(dec, member, reversed, pArc));
   dec->arcs[*pArc].element = isTree ? MARKER_ROW_ELEMENT : MARKER_COLUMN_ELEMENT;
   dec->arcs[*pArc].childMember = child;

   addArcToMemberArcList(dec, *pArc, member);
   return SCIP_OKAY;
}

static SCIP_RETCODE createParentMarker(
   SCIP_NETMATDECDATA* dec,
   spqr_member member,
   SCIP_Bool isTree,
   spqr_member parent,
   spqr_arc parentMarker,
   spqr_arc* arc,
   SCIP_Bool reversed
)
{

   SCIP_CALL(createArc(dec, member, reversed, arc));
   dec->arcs[*arc].element = isTree ? MARKER_ROW_ELEMENT : MARKER_COLUMN_ELEMENT;

   addArcToMemberArcList(dec, *arc, member);

   dec->members[member].parentMember = parent;
   dec->members[member].markerOfParent = parentMarker;
   dec->members[member].markerToParent = *arc;
   return SCIP_OKAY;
}

static SCIP_RETCODE createMarkerPair(
   SCIP_NETMATDECDATA* dec,
   spqr_member parentMember,
   spqr_member childMember,
   SCIP_Bool parentIsTree,
   SCIP_Bool childMarkerReversed,
   SCIP_Bool parentMarkerReversed
)
{
   spqr_arc parentToChildMarker = SPQR_INVALID_ARC;
   SCIP_CALL(
      createChildMarker(dec, parentMember, childMember, parentIsTree, &parentToChildMarker, childMarkerReversed));

   spqr_arc childToParentMarker = SPQR_INVALID_ARC;
   SCIP_CALL(
      createParentMarker(dec, childMember, !parentIsTree, parentMember, parentToChildMarker, &childToParentMarker,
                         parentMarkerReversed));

   return SCIP_OKAY;
}

static SCIP_RETCODE createMarkerPairWithReferences(
   SCIP_NETMATDECDATA* dec,
   spqr_member parentMember,
   spqr_member childMember,
   SCIP_Bool parentIsTree,
   SCIP_Bool childMarkerReversed,
   SCIP_Bool parentMarkerReversed,
   spqr_arc* parentToChild,
   spqr_arc* childToParent
)
{
   SCIP_CALL(createChildMarker(dec, parentMember, childMember, parentIsTree, parentToChild, childMarkerReversed));
   SCIP_CALL(createParentMarker(dec, childMember, !parentIsTree, parentMember, *parentToChild, childToParent,
                                parentMarkerReversed));

   return SCIP_OKAY;
}

static void moveArcToNewMember(
   SCIP_NETMATDECDATA* dec,
   spqr_arc arc,
   spqr_member oldMember,
   spqr_member newMember
)
{
   assert(SPQRarcIsValid(arc));
   assert(arc < dec->memArcs);
   assert(dec);

   assert(memberIsRepresentative(dec, oldMember));
   assert(memberIsRepresentative(dec, newMember));
   //Need to change the arc's member, remove it from the current member list and add it to the new member list
   assert(findArcMemberNoCompression(dec, arc) == oldMember);

   removeArcFromMemberArcList(dec, arc, oldMember);
   addArcToMemberArcList(dec, arc, newMember);

   dec->arcs[arc].member = newMember;

   //If this arc has a childMember, update the information correctly!
   spqr_member childMember = dec->arcs[arc].childMember;
   if( SPQRmemberIsValid(childMember))
   {
      spqr_member childRepresentative = findArcChildMember(dec, arc);
      dec->members[childRepresentative].parentMember = newMember;
   }
   //If this arc is a marker to the parent, update the child arc marker of the parent to reflect the move
   if( dec->members[oldMember].markerToParent == arc )
   {
      dec->members[newMember].markerToParent = arc;
      dec->members[newMember].parentMember = dec->members[oldMember].parentMember;
      dec->members[newMember].markerOfParent = dec->members[oldMember].markerOfParent;

      assert(findArcChildMemberNoCompression(dec, dec->members[oldMember].markerOfParent) == oldMember);
      dec->arcs[dec->members[oldMember].markerOfParent].childMember = newMember;
   }
}

static void mergeMemberArcList(
   SCIP_NETMATDECDATA* dec,
   spqr_member toMergeInto,
   spqr_member toRemove
)
{
   spqr_arc firstIntoArc = getFirstMemberArc(dec, toMergeInto);
   spqr_arc firstFromArc = getFirstMemberArc(dec, toRemove);
   assert(SPQRarcIsValid(firstIntoArc));
   assert(SPQRarcIsValid(firstFromArc));

   spqr_arc lastIntoArc = getPreviousMemberArc(dec, firstIntoArc);
   spqr_arc lastFromArc = getPreviousMemberArc(dec, firstFromArc);

   //Relink linked lists to merge them effectively
   dec->arcs[firstIntoArc].arcListNode.previous = lastFromArc;
   dec->arcs[lastIntoArc].arcListNode.next = firstFromArc;
   dec->arcs[firstFromArc].arcListNode.previous = lastIntoArc;
   dec->arcs[lastFromArc].arcListNode.next = firstIntoArc;

   //Clean up old
   dec->members[toMergeInto].numArcs += dec->members[toRemove].numArcs;
   dec->members[toRemove].numArcs = 0;
   dec->members[toRemove].firstArc = SPQR_INVALID_ARC;
}

static void changeLoopToSeries(
   SCIP_NETMATDECDATA* dec,
   spqr_member member
)
{
   assert(SPQRmemberIsValid(member));
   assert(member < dec->memMembers);
   assert(dec);
   assert(( getMemberType(dec, member) == SPQR_MEMBERTYPE_PARALLEL ||
            getMemberType(dec, member) == SPQR_MEMBERTYPE_SERIES ||
            getMemberType(dec, member) == SPQR_MEMBERTYPE_LOOP ) &&
          getNumMemberArcs(dec, member) == 2);
   assert(memberIsRepresentative(dec, member));
   dec->members[member].type = SPQR_MEMBERTYPE_SERIES;
}

static void changeLoopToParallel(
   SCIP_NETMATDECDATA* dec,
   spqr_member member
)
{
   assert(SPQRmemberIsValid(member));
   assert(member < dec->memMembers);
   assert(dec);
   assert(( getMemberType(dec, member) == SPQR_MEMBERTYPE_PARALLEL ||
            getMemberType(dec, member) == SPQR_MEMBERTYPE_SERIES ||
            getMemberType(dec, member) == SPQR_MEMBERTYPE_LOOP ) &&
          getNumMemberArcs(dec, member) == 2);
   assert(memberIsRepresentative(dec, member));
   dec->members[member].type = SPQR_MEMBERTYPE_PARALLEL;
}

static SCIP_Bool netMatDecDataIsMinimal(const SCIP_NETMATDECDATA* dec)
{
   //Relies on parents/children etc. being set correctly in the tree
   SCIP_Bool isMinimal = TRUE;
   for( spqr_member member = 0; member < dec->numMembers; ++member )
   {
      if( !memberIsRepresentative(dec, member) || getMemberType(dec, member) == SPQR_MEMBERTYPE_UNASSIGNED )
      {
         continue;
      }
      spqr_member memberParent = findMemberParentNoCompression(dec, member);
      if( SPQRmemberIsValid(memberParent))
      {
         assert(memberIsRepresentative(dec, memberParent));
         SPQRMemberType memberType = getMemberType(dec, member);
         SPQRMemberType parentType = getMemberType(dec, memberParent);
         if( memberType == parentType && memberType != SPQR_MEMBERTYPE_RIGID )
         {
            isMinimal = FALSE;
            break;
         }
      }
   }
   return isMinimal;
}

static void decreaseNumConnectedComponents(
   SCIP_NETMATDECDATA* dec,
   int by
)
{
   dec->numConnectedComponents -= by;
   assert(dec->numConnectedComponents >= 1);
}

static void reorderComponent(
   SCIP_NETMATDECDATA* dec,
   spqr_member newRoot
)
{
   assert(dec);
   assert(memberIsRepresentative(dec, newRoot));
   //If the newRoot has no parent, it is already the root, so then there's no need to reorder.
   if( SPQRmemberIsValid(dec->members[newRoot].parentMember))
   {
      spqr_member member = findMemberParent(dec, newRoot);
      spqr_member newParent = newRoot;
      spqr_arc newMarkerToParent = dec->members[newRoot].markerOfParent;
      spqr_arc markerOfNewParent = dec->members[newRoot].markerToParent;

      //Recursively update the parent
      do
      {
         assert(SPQRmemberIsValid(member));
         assert(SPQRmemberIsValid(newParent));
         spqr_member oldParent = findMemberParent(dec, member);
         spqr_arc oldMarkerToParent = dec->members[member].markerToParent;
         spqr_arc oldMarkerOfParent = dec->members[member].markerOfParent;

         dec->members[member].markerToParent = newMarkerToParent;
         dec->members[member].markerOfParent = markerOfNewParent;
         dec->members[member].parentMember = newParent;
         dec->arcs[markerOfNewParent].childMember = member;
         dec->arcs[newMarkerToParent].childMember = SPQR_INVALID_MEMBER;

         if( SPQRmemberIsValid(oldParent))
         {
            newParent = member;
            member = oldParent;
            newMarkerToParent = oldMarkerOfParent;
            markerOfNewParent = oldMarkerToParent;
         } else
         {
            break;
         }
      } while( TRUE ); /*lint !e506*/
      dec->members[newRoot].parentMember = SPQR_INVALID_MEMBER;
      dec->members[newRoot].markerToParent = SPQR_INVALID_ARC;
      dec->members[newRoot].markerOfParent = SPQR_INVALID_ARC;
   }
}


static void netMatDecDataRemoveComponent(
   SCIP_NETMATDECDATA* dec,
   const int* componentRows,
   int numRows,
   const int* componentCols,
   int numCols
)
{
   //The below just removes the 'link' but not the internal datastructures.
   //This is sufficient for our purposes, as long as we do not re-introduce any of the 'negated' rows/columns back into the decomposition.

   for( int i = 0; i < numRows; ++i )
   {
      spqr_row row = componentRows[i];
      if( SPQRarcIsValid(dec->rowArcs[row]))
      {
         dec->rowArcs[row] = SPQR_INVALID_ARC;
      }
   }

   for( int i = 0; i < numCols; ++i )
   {
      spqr_col col = componentCols[i];
      if( SPQRarcIsValid(dec->columnArcs[col]))
      {
         dec->columnArcs[col] = SPQR_INVALID_ARC;
      }
   }
}

#ifdef SCIP_DEBUG
                                                                                                                        //Debugging functions to print the decomposition
static char typeToChar(SPQRMemberType type)
{
   switch (type)
   {
      case SPQR_MEMBERTYPE_RIGID:
         return 'R';
      case SPQR_MEMBERTYPE_PARALLEL:
         return 'P';
      case SPQR_MEMBERTYPE_SERIES:
         return 'S';
      case SPQR_MEMBERTYPE_LOOP:
         return 'L';
      default:
         return '?';
   }
}

static void arcToDot(FILE *stream, const SCIP_NETMATDEC *dec,
                     spqr_arc arc, unsigned long dot_head, unsigned long dot_tail, SCIP_Bool useElementNames)
{
   assert(SPQRarcIsValid(arc));
   spqr_member member = findArcMemberNoCompression(dec, arc);
   SPQRMemberType member_type = getMemberType(dec, member);
   char type = typeToChar(member_type);
   const char *color = arcIsTree(dec, arc) ? ",color=red" : ",color=blue";

   int arc_name = arc;

   if (markerToParent(dec, member) == arc)
   {
      if (useElementNames)
      {
         arc_name = -1;
      }
      fprintf(stream, "    %c_%d_%lu -> %c_p_%d [label=\"%d\",style=dashed%s];\n", type, member, dot_tail, type, member, arc_name, color);
      fprintf(stream, "    %c_p_%d -> %c_%d_%lu [label=\"%d\",style=dashed%s];\n", type, member, type, member, dot_head, arc_name, color);
      fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_tail);
      fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_head);
      fprintf(stream, "    %c_p_%d [style=dashed];\n", type, member);
   } else if (arcIsMarker(dec, arc))
   {
      spqr_member child = findArcChildMemberNoCompression(dec, arc);
      char childType = typeToChar(getMemberType(dec, child));
      if (useElementNames)
      {
         arc_name = -1;
      }
      fprintf(stream, "    %c_%d_%lu -> %c_c_%d [label=\"%d\",style=dotted%s];\n", type, member, dot_tail, type, child, arc_name, color);
      fprintf(stream, "    %c_c_%d -> %c_%d_%lu [label=\"%d\",style=dotted%s];\n", type, child, type, member, dot_head, arc_name, color);
      fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_tail);
      fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_head);
      fprintf(stream, "    %c_c_%d [style=dotted];\n", type, child);
      fprintf(stream, "    %c_p_%d -> %c_c_%d [style=dashed,dir=forward];\n", childType, child, type, child);
   } else
   {
      if (useElementNames)
      {
         spqr_element element = dec->arcs[arc].element;
         if (SPQRelementIsRow(element))
         {
            arc_name = (int) SPQRelementToRow(element);
         } else
         {
            arc_name = (int) SPQRelementToColumn(element);
         }
      }

      fprintf(stream, "    %c_%d_%lu -> %c_%d_%lu [label=\"%d \",style=bold%s];\n", type, member, dot_tail, type, member, dot_head,
              arc_name, color);
      fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_tail);
      fprintf(stream, "    %c_%d_%lu [shape=box];\n", type, member, dot_head);
   }
}

static void decompositionToDot(FILE *stream, const SCIP_NETMATDEC *dec, SCIP_Bool useElementNames)
{
   fprintf(stream, "//decomposition\ndigraph decomposition{\n   compound = TRUE;\n");
   for (spqr_member member = 0; member < dec->numMembers; ++member)
   {
      if (!memberIsRepresentative(dec, member)) continue;
      fprintf(stream, "   subgraph member_%d{\n", member);
      switch (getMemberType(dec, member))
      {
         case SPQR_MEMBERTYPE_RIGID: {
            spqr_arc first_arc = getFirstMemberArc(dec, member);
            spqr_arc arc = first_arc;
            do
            {
               unsigned long arcHead = (unsigned long) findEffectiveArcHeadNoCompression(dec, arc);
               unsigned long arcTail = (unsigned long) findEffectiveArcTailNoCompression(dec, arc);
               arcToDot(stream, dec, arc, arcHead, arcTail, useElementNames);
               arc = getNextMemberArc(dec, arc);
            } while (arc != first_arc);
            break;
         }
         case SPQR_MEMBERTYPE_LOOP:
         case SPQR_MEMBERTYPE_PARALLEL: {
            spqr_arc first_arc = getFirstMemberArc(dec, member);
            spqr_arc arc = first_arc;
            do
            {
               if (arcIsReversedNonRigid(dec, arc))
               {
                  arcToDot(stream, dec, arc, 1, 0, useElementNames);
               } else
               {
                  arcToDot(stream, dec, arc, 0, 1, useElementNames);
               }
               arc = getNextMemberArc(dec, arc);
            } while (arc != first_arc);
            break;
         }
         case SPQR_MEMBERTYPE_SERIES: {
            unsigned long i = 0;
            unsigned long num_member_arcs = (unsigned long) getNumMemberArcs(dec, member);
            spqr_arc first_arc = getFirstMemberArc(dec, member);
            spqr_arc arc = first_arc;
            do {
               unsigned long head = i;
               unsigned long tail = (i + 1) % num_member_arcs;
               if (arcIsReversedNonRigid(dec, arc))
               {
                  unsigned long temp = head;
                  head = tail;
                  tail = temp;
               }
               arcToDot(stream, dec, arc, head, tail, useElementNames);
               arc = getNextMemberArc(dec, arc);
               i++;
            } while (arc != first_arc);
            break;
         }
         case SPQR_MEMBERTYPE_UNASSIGNED:
            break;
      }
      fprintf(stream, "   }\n");
   }
   fprintf(stream, "}\n");
}
#endif

static SCIP_RETCODE mergeGivenMemberIntoParent(
   SCIP_NETMATDECDATA* dec,
   spqr_member member,
   spqr_member parent,
   spqr_arc parentToChild,
   spqr_arc childToParent,
   SCIP_Bool headToHead,
   spqr_member* mergedMember
)
{
   assert(dec);
   assert(SPQRmemberIsValid(member));
   assert(memberIsRepresentative(dec, member));
   assert(SPQRmemberIsValid(parent));
   assert(memberIsRepresentative(dec, parent));
   assert(findMemberParentNoCompression(dec, member) == parent);
   assert(markerOfParent(dec, member) == parentToChild);
   assert(markerToParent(dec, member) == childToParent);

   removeArcFromMemberArcList(dec, parentToChild, parent);
   removeArcFromMemberArcList(dec, childToParent, member);

   spqr_node parentTail = findEffectiveArcTail(dec, parentToChild) ;
   spqr_node parentHead = findEffectiveArcHead(dec, parentToChild);
   spqr_node childTail =  findEffectiveArcTail(dec, childToParent);
   spqr_node childHead = findEffectiveArcHead(dec, childToParent);
   spqr_node parentArcNodes[2] = { parentTail, parentHead };
   spqr_node childArcNodes[2] = { childTail, childHead };

   clearArcHeadAndTail(dec, parentToChild);
   clearArcHeadAndTail(dec, childToParent);

   spqr_node first = childArcNodes[headToHead ? 0 : 1];
   spqr_node second = childArcNodes[headToHead ? 1 : 0];
   {
      spqr_node newNode = mergeNodes(dec, parentArcNodes[0], first);
      spqr_node toRemoveFrom = newNode == first ? parentArcNodes[0] : first;
      mergeNodeArcList(dec, newNode, toRemoveFrom);
   }
   {
      spqr_node newNode = mergeNodes(dec, parentArcNodes[1], second);
      spqr_node toRemoveFrom = newNode == second ? parentArcNodes[1] : second;
      mergeNodeArcList(dec, newNode, toRemoveFrom);
   }


   spqr_member newMember = mergeMembers(dec, member, parent);
   spqr_member toRemoveFrom = newMember == member ? parent : member;
   mergeMemberArcList(dec, newMember, toRemoveFrom);
   if( toRemoveFrom == parent )
   {
      updateMemberParentInformation(dec, newMember, toRemoveFrom);
   }
   updateMemberType(dec, newMember, SPQR_MEMBERTYPE_RIGID);
   *mergedMember = newMember;
   return SCIP_OKAY;
}

static int maxValue(
   int a,
   int b
)
{
   return ( a > b ) ? a : b;
}

typedef int path_arc_id;
#define INVALID_PATH_ARC (-1)

static SCIP_Bool pathArcIsInvalid(const path_arc_id arc)
{
   return arc < 0;
}

static SCIP_Bool pathArcIsValid(const path_arc_id arc)
{
   return !pathArcIsInvalid(arc);
}

typedef struct
{
   spqr_arc arc;
   spqr_node arcHead;//These can be used in various places to prevent additional find()'s
   spqr_node arcTail;
   path_arc_id nextMember;
   path_arc_id nextOverall;
   SCIP_Bool reversed;
} PathArcListNode;

typedef int reduced_member_id;
#define INVALID_REDUCED_MEMBER (-1)

static SCIP_Bool reducedMemberIsInvalid(const reduced_member_id id)
{
   return id < 0;
}

static SCIP_Bool reducedMemberIsValid(const reduced_member_id id)
{
   return !reducedMemberIsInvalid(id);
}

typedef int children_idx;

typedef enum
{
   REDUCEDMEMBER_TYPE_UNASSIGNED = 0,
   REDUCEDMEMBER_TYPE_CYCLE = 1,
   REDUCEDMEMBER_TYPE_MERGED = 2,
   REDUCEDMEMBER_TYPE_NOT_NETWORK = 3
} ReducedMemberType;

typedef enum
{
   INTO_HEAD = 0,
   INTO_TAIL = 1,
   OUT_HEAD = 2,
   OUT_TAIL = 3
} MemberPathType;

static SCIP_Bool isInto(MemberPathType type)
{
   return type == INTO_HEAD || type == INTO_TAIL;
}

static SCIP_Bool isHead(MemberPathType type)
{
   return type == INTO_HEAD || type == OUT_HEAD;
}

typedef struct
{
   spqr_member member;
   spqr_member rootMember;
   int depth;
   ReducedMemberType type;
   reduced_member_id parent;

   children_idx firstChild;
   children_idx numChildren;

   path_arc_id firstPathArc;
   int numPathArcs;

   SCIP_Bool reverseArcs;
   spqr_node rigidPathStart;
   spqr_node rigidPathEnd;

   SCIP_Bool pathBackwards;

   int numPropagatedChildren;
   int componentIndex;

   MemberPathType pathType;
   reduced_member_id nextPathMember;
   SCIP_Bool nextPathMemberIsParent;
   spqr_arc pathSourceArc;
   spqr_arc pathTargetArc;
} SPQRColReducedMember;

typedef struct
{
   int rootDepth;
   reduced_member_id root;

   reduced_member_id pathEndMembers[2];
   int numPathEndMembers;
} SPQRColReducedComponent;

typedef struct
{
   reduced_member_id reducedMember;
   reduced_member_id rootDepthMinimizer;
} MemberInfo;

typedef struct
{
   spqr_member member;
} CreateReducedMembersCallstack;

typedef struct
{
   SCIP_Bool remainsNetwork;

   SPQRColReducedMember* reducedMembers;
   int memReducedMembers;
   int numReducedMembers;

   SPQRColReducedComponent* reducedComponents;
   int memReducedComponents;
   int numReducedComponents;

   MemberInfo* memberInformation;
   int memMemberInformation;
   int numMemberInformation;

   reduced_member_id* childrenStorage;
   int memChildrenStorage;
   int numChildrenStorage;

   PathArcListNode* pathArcs;
   int memPathArcs;
   int numPathArcs;
   path_arc_id firstOverallPathArc;

   int* nodeInPathDegree;
   int* nodeOutPathDegree;
   int memNodePathDegree;

   SCIP_Bool* arcInPath;
   SCIP_Bool* arcInPathReversed;
   int memArcsInPath;

   CreateReducedMembersCallstack* createReducedMembersCallStack;
   int memCreateReducedMembersCallStack;

   spqr_col newColIndex;

   spqr_row* newRowArcs;
   SCIP_Bool* newRowArcReversed;
   int memNewRowArcs;
   int numNewRowArcs;

   spqr_arc* decompositionRowArcs;
   SCIP_Bool* decompositionArcReversed;
   int memDecompositionRowArcs;
   int numDecompositionRowArcs;

   spqr_member* leafMembers;
   int numLeafMembers;
   int memLeafMembers;
} SCIP_NETCOLADD;

static void cleanupPreviousIteration(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol
)
{
   assert(dec);
   assert(newCol);

   path_arc_id pathArc = newCol->firstOverallPathArc;
   while( pathArcIsValid(pathArc))
   {
      spqr_node head = newCol->pathArcs[pathArc].arcHead;
      spqr_node tail = newCol->pathArcs[pathArc].arcTail;
      if( SPQRnodeIsValid(head))
      {
         newCol->nodeInPathDegree[head] = 0;
      }
      if( SPQRnodeIsValid(tail))
      {
         newCol->nodeOutPathDegree[tail] = 0;
      }

      spqr_arc arc = newCol->pathArcs[pathArc].arc;
      if( arc < newCol->memArcsInPath )
      {
         newCol->arcInPath[arc] = FALSE;
         newCol->arcInPathReversed[arc] = FALSE;
      }
      pathArc = newCol->pathArcs[pathArc].nextOverall;
   }
#ifndef NDEBUG
   for( int i = 0; i < newCol->memArcsInPath; ++i )
   {
      assert(newCol->arcInPath[i] == FALSE);
      assert(newCol->arcInPathReversed[i] == FALSE);
   }

   for( int i = 0; i < newCol->memNodePathDegree; ++i )
   {
      assert(newCol->nodeInPathDegree[i] == 0);
      assert(newCol->nodeOutPathDegree[i] == 0);
   }
#endif

   newCol->firstOverallPathArc = INVALID_PATH_ARC;
   newCol->numPathArcs = 0;
}

static SCIP_RETCODE SCIPnetcoladdCreate(
   SCIP* scip,
   SCIP_NETCOLADD** pcoladd
)
{
   assert(scip);

   SCIP_CALL(SCIPallocBlockMemory(scip, pcoladd));
   SCIP_NETCOLADD* newCol = *pcoladd;

   newCol->remainsNetwork = FALSE;
   newCol->reducedMembers = NULL;
   newCol->memReducedMembers = 0;
   newCol->numReducedMembers = 0;

   newCol->reducedComponents = NULL;
   newCol->memReducedComponents = 0;
   newCol->numReducedComponents = 0;

   newCol->memberInformation = NULL;
   newCol->memMemberInformation = 0;
   newCol->numMemberInformation = 0;

   newCol->childrenStorage = NULL;
   newCol->memChildrenStorage = 0;
   newCol->numChildrenStorage = 0;

   newCol->pathArcs = NULL;
   newCol->memPathArcs = 0;
   newCol->numPathArcs = 0;
   newCol->firstOverallPathArc = INVALID_PATH_ARC;

   newCol->nodeInPathDegree = NULL;
   newCol->nodeOutPathDegree = NULL;
   newCol->memNodePathDegree = 0;

   newCol->arcInPath = NULL;
   newCol->arcInPathReversed = NULL;
   newCol->memArcsInPath = 0;

   newCol->createReducedMembersCallStack = NULL;
   newCol->memCreateReducedMembersCallStack = 0;

   newCol->newColIndex = SPQR_INVALID_COL;

   newCol->newRowArcs = NULL;
   newCol->newRowArcReversed = NULL;
   newCol->memNewRowArcs = 0;
   newCol->numNewRowArcs = 0;

   newCol->decompositionRowArcs = NULL;
   newCol->decompositionArcReversed = NULL;
   newCol->memDecompositionRowArcs = 0;
   newCol->numDecompositionRowArcs = 0;

   newCol->leafMembers = NULL;
   newCol->numLeafMembers = 0;
   newCol->memLeafMembers = 0;

   return SCIP_OKAY;
}

static void SCIPnetcoladdFree(
   SCIP* scip,
   SCIP_NETCOLADD** pcoladd
)
{
   assert(scip);
   SCIP_NETCOLADD* newCol = *pcoladd;
   SCIPfreeBlockMemoryArray(scip, &newCol->decompositionRowArcs, newCol->memDecompositionRowArcs);
   SCIPfreeBlockMemoryArray(scip, &newCol->decompositionArcReversed, newCol->memDecompositionRowArcs);
   SCIPfreeBlockMemoryArray(scip, &newCol->newRowArcs, newCol->memNewRowArcs);
   SCIPfreeBlockMemoryArray(scip, &newCol->newRowArcReversed, newCol->memNewRowArcs);
   SCIPfreeBlockMemoryArray(scip, &newCol->createReducedMembersCallStack, newCol->memCreateReducedMembersCallStack);
   SCIPfreeBlockMemoryArray(scip, &newCol->arcInPath, newCol->memArcsInPath);
   SCIPfreeBlockMemoryArray(scip, &newCol->arcInPathReversed, newCol->memArcsInPath);
   SCIPfreeBlockMemoryArray(scip, &newCol->nodeInPathDegree, newCol->memNodePathDegree);
   SCIPfreeBlockMemoryArray(scip, &newCol->nodeOutPathDegree, newCol->memNodePathDegree);
   SCIPfreeBlockMemoryArray(scip, &newCol->pathArcs, newCol->memPathArcs);
   SCIPfreeBlockMemoryArray(scip, &newCol->childrenStorage, newCol->memChildrenStorage);
   SCIPfreeBlockMemoryArray(scip, &newCol->memberInformation, newCol->memMemberInformation);
   SCIPfreeBlockMemoryArray(scip, &newCol->reducedComponents, newCol->memReducedComponents);
   SCIPfreeBlockMemoryArray(scip, &newCol->reducedMembers, newCol->memReducedMembers);
   SCIPfreeBlockMemoryArray(scip, &newCol->leafMembers, newCol->memLeafMembers);

   SCIPfreeBlockMemory(scip, pcoladd);
}


static reduced_member_id createReducedMembersToRoot(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   const spqr_member firstMember
)
{
   assert(SPQRmemberIsValid(firstMember));

   CreateReducedMembersCallstack* callstack = newCol->createReducedMembersCallStack;
   callstack[0].member = firstMember;
   int callDepth = 0;

   while( callDepth >= 0 )
   {
      spqr_member member = callstack[callDepth].member;
      reduced_member_id reducedMember = newCol->memberInformation[member].reducedMember;

      SCIP_Bool reducedValid = reducedMemberIsValid(reducedMember);
      if( !reducedValid )
      {
         //reduced member was not yet created; we create it
         reducedMember = newCol->numReducedMembers;

         SPQRColReducedMember* reducedMemberData = &newCol->reducedMembers[reducedMember];
         ++newCol->numReducedMembers;

         reducedMemberData->member = member;
         reducedMemberData->numChildren = 0;

         reducedMemberData->type = REDUCEDMEMBER_TYPE_UNASSIGNED;
         reducedMemberData->numPropagatedChildren = 0;
         reducedMemberData->firstPathArc = INVALID_PATH_ARC;
         reducedMemberData->numPathArcs = 0;
         reducedMemberData->rigidPathStart = SPQR_INVALID_NODE;
         reducedMemberData->rigidPathEnd = SPQR_INVALID_NODE;

         reducedMemberData->componentIndex = -1;
         //The children are set later

         newCol->memberInformation[member].reducedMember = reducedMember;
         assert(memberIsRepresentative(dec, member));
         spqr_member parentMember = findMemberParent(dec, member);

         if( SPQRmemberIsValid(parentMember))
         {
            //recursive call to parent member
            ++callDepth;
            assert(callDepth < newCol->memCreateReducedMembersCallStack);
            callstack[callDepth].member = parentMember;
            continue;

         } else
         {
            //we found a new reduced decomposition component

            reducedMemberData->parent = INVALID_REDUCED_MEMBER;
            reducedMemberData->depth = 0;
            reducedMemberData->rootMember = member;
            reducedMemberData->componentIndex = newCol->numReducedComponents;

            assert(newCol->numReducedComponents < newCol->memReducedComponents);
            newCol->reducedComponents[newCol->numReducedComponents].root = reducedMember;
            newCol->reducedComponents[newCol->numReducedComponents].numPathEndMembers = 0;
            newCol->reducedComponents[newCol->numReducedComponents].pathEndMembers[0] = INVALID_REDUCED_MEMBER;
            newCol->reducedComponents[newCol->numReducedComponents].pathEndMembers[1] = INVALID_REDUCED_MEMBER;
            ++newCol->numReducedComponents;
         }
      }
      if( reducedValid )
      {
         assert(reducedMember < newCol->numReducedMembers);
         //Reduced member was already created in earlier call
         //update the depth of the root if appropriate
         reduced_member_id* depthMinimizer = &newCol->memberInformation[newCol->reducedMembers[reducedMember].rootMember].rootDepthMinimizer;
         if( reducedMemberIsInvalid(*depthMinimizer) ||
             newCol->reducedMembers[reducedMember].depth < newCol->reducedMembers[*depthMinimizer].depth )
         {
            *depthMinimizer = reducedMember;
         }
      }
      while( TRUE ) /*lint !e716*/
      {
         --callDepth;
         if( callDepth < 0 ) break;
         spqr_member parentMember = callstack[callDepth + 1].member;
         reduced_member_id parentReducedMember = newCol->memberInformation[parentMember].reducedMember;
         spqr_member currentMember = callstack[callDepth].member;
         reduced_member_id currentReducedMember = newCol->memberInformation[currentMember].reducedMember;

         SPQRColReducedMember* parentReducedMemberData = &newCol->reducedMembers[parentReducedMember];
         SPQRColReducedMember* reducedMemberData = &newCol->reducedMembers[currentReducedMember];

         reducedMemberData->parent = parentReducedMember;
         reducedMemberData->depth = parentReducedMemberData->depth + 1;
         reducedMemberData->rootMember = parentReducedMemberData->rootMember;
         //ensure that all newly created reduced members are pointing to the correct component
         assert(parentReducedMemberData->componentIndex >= 0);
         reducedMemberData->componentIndex = parentReducedMemberData->componentIndex;

         newCol->reducedMembers[parentReducedMember].numChildren++;
      }
   }

   reduced_member_id returnedMember = newCol->memberInformation[callstack[0].member].reducedMember;
   return returnedMember;
}

static SCIP_RETCODE constructReducedDecomposition(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol
)
{
   assert(dec);
   assert(newCol);
#ifndef NDEBUG
   for( int i = 0; i < newCol->memMemberInformation; ++i )
   {
      assert(reducedMemberIsInvalid(newCol->memberInformation[i].reducedMember));
   }
#endif
   newCol->numReducedComponents = 0;
   newCol->numReducedMembers = 0;
   if( newCol->numDecompositionRowArcs == 0 )
   {//Early return in case the reduced decomposition will be empty
      return SCIP_OKAY;
   }
   assert(newCol->numReducedMembers == 0);
   assert(newCol->numReducedComponents == 0);

   int newSize = largestMemberID(dec);//Is this sufficient?
   if( newSize > newCol->memReducedMembers )
   {
      int newArraySize = maxValue(2 * newCol->memReducedMembers, newSize);
      SCIP_CALL(
         SCIPreallocBlockMemoryArray(dec->env, &newCol->reducedMembers, newCol->memReducedMembers, newArraySize));
      newCol->memReducedMembers = newArraySize;
   }
   if( newSize > newCol->memMemberInformation )
   {
      int updatedSize = maxValue(2 * newCol->memMemberInformation, newSize);
      SCIP_CALL(
         SCIPreallocBlockMemoryArray(dec->env, &newCol->memberInformation, newCol->memMemberInformation, updatedSize));
      for( int i = newCol->memMemberInformation; i < updatedSize; ++i )
      {
         newCol->memberInformation[i].reducedMember = INVALID_REDUCED_MEMBER;
         newCol->memberInformation[i].rootDepthMinimizer = INVALID_REDUCED_MEMBER;
      }
      newCol->memMemberInformation = updatedSize;
   }

   int numComponents = numConnectedComponents(dec);
   if( numComponents > newCol->memReducedComponents )
   {
      int updatedSize = maxValue(2 * newCol->memReducedComponents, numComponents);
      SCIP_CALL(
         SCIPreallocBlockMemoryArray(dec->env, &newCol->reducedComponents, newCol->memReducedComponents, updatedSize));
      newCol->memReducedComponents = updatedSize;
   }

   int numMembers = getNumMembers(dec);
   if( newCol->memCreateReducedMembersCallStack < numMembers )
   {
      int updatedSize = maxValue(2 * newCol->memCreateReducedMembersCallStack, numMembers);
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newCol->createReducedMembersCallStack,
                                            newCol->memCreateReducedMembersCallStack, updatedSize));
      newCol->memCreateReducedMembersCallStack = updatedSize;
   }

   //Create the reduced members (recursively)
   for( int i = 0; i < newCol->numDecompositionRowArcs; ++i )
   {
      assert(i < newCol->memDecompositionRowArcs);
      spqr_arc arc = newCol->decompositionRowArcs[i];
      spqr_member arcMember = findArcMember(dec, arc);
      reduced_member_id reducedMember = createReducedMembersToRoot(dec, newCol, arcMember);
      reduced_member_id* depthMinimizer = &newCol->memberInformation[newCol->reducedMembers[reducedMember].rootMember].rootDepthMinimizer;
      if( reducedMemberIsInvalid(*depthMinimizer))
      {
         *depthMinimizer = reducedMember;
      }
   }

   //Set the reduced roots according to the root depth minimizers
   for( int i = 0; i < newCol->numReducedComponents; ++i )
   {
      SPQRColReducedComponent* component = &newCol->reducedComponents[i];
      spqr_member rootMember = newCol->reducedMembers[component->root].member;
      reduced_member_id reducedMinimizer = newCol->memberInformation[rootMember].rootDepthMinimizer;
      component->rootDepth = newCol->reducedMembers[reducedMinimizer].depth;
      component->root = reducedMinimizer;

      //This simplifies code further down which does not need to be component-aware; just pretend that the reduced member is the new root.
      newCol->reducedMembers[component->root].parent = INVALID_REDUCED_MEMBER;
      assert(memberIsRepresentative(dec, rootMember));
   }

   //update the children array
   int numTotalChildren = 0;
   for( int i = 0; i < newCol->numReducedMembers; ++i )
   {
      SPQRColReducedMember* reducedMember = &newCol->reducedMembers[i];
      reduced_member_id minimizer = newCol->memberInformation[reducedMember->rootMember].rootDepthMinimizer;
      if( reducedMember->depth >= newCol->reducedMembers[minimizer].depth )
      {
         reducedMember->firstChild = numTotalChildren;
         numTotalChildren += reducedMember->numChildren;
         reducedMember->numChildren = 0;
      }
   }

   if( newCol->memChildrenStorage < numTotalChildren )
   {
      int newMemSize = maxValue(newCol->memChildrenStorage * 2, numTotalChildren);
      SCIP_CALL(
         SCIPreallocBlockMemoryArray(dec->env, &newCol->childrenStorage, newCol->memChildrenStorage, newMemSize));
      newCol->memChildrenStorage = newMemSize;
   }
   newCol->numChildrenStorage = numTotalChildren;

   //Fill up the children array`
   for( reduced_member_id reducedMember = 0; reducedMember < newCol->numReducedMembers; ++reducedMember )
   {
      SPQRColReducedMember* reducedMemberData = &newCol->reducedMembers[reducedMember];
      if( reducedMemberData->depth <=
          newCol->reducedMembers[newCol->memberInformation[reducedMemberData->rootMember].rootDepthMinimizer].depth )
      {
         continue;
      }
      spqr_member parentMember = findMemberParent(dec, reducedMemberData->member);
      reduced_member_id parentReducedMember = SPQRmemberIsValid(parentMember)
                                              ? newCol->memberInformation[parentMember].reducedMember
                                              : INVALID_REDUCED_MEMBER;
      if( reducedMemberIsValid(parentReducedMember))
      {
         SPQRColReducedMember* parentReducedMemberData = &newCol->reducedMembers[parentReducedMember];
         newCol->childrenStorage[parentReducedMemberData->firstChild +
                                 parentReducedMemberData->numChildren] = reducedMember;
         ++parentReducedMemberData->numChildren;
      }
   }

   //Clean up the root depth minimizers.
   for( int i = 0; i < newCol->numReducedMembers; ++i )
   {
      SPQRColReducedMember* reducedMember = &newCol->reducedMembers[i];
      assert(reducedMember);
      spqr_member rootMember = reducedMember->rootMember;
      assert(rootMember >= 0);
      assert(rootMember < dec->memMembers);
      newCol->memberInformation[rootMember].rootDepthMinimizer = INVALID_REDUCED_MEMBER;
   }

   return SCIP_OKAY;
}

static void cleanUpMemberInformation(SCIP_NETCOLADD* newCol)
{
   //This loop is at the end as memberInformation is also used to assign the cut arcs during propagation
   //Clean up the memberInformation array
   for( int i = 0; i < newCol->numReducedMembers; ++i )
   {
      newCol->memberInformation[newCol->reducedMembers[i].member].reducedMember = INVALID_REDUCED_MEMBER;
   }
#ifndef NDEBUG
   for( int i = 0; i < newCol->memMemberInformation; ++i )
   {
      assert(reducedMemberIsInvalid(newCol->memberInformation[i].reducedMember));
   }
#endif
}

static void createPathArc(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   const spqr_arc arc,
   const reduced_member_id reducedMember,
   SCIP_Bool reversed
)
{
   assert(dec);
   assert(newCol);

   path_arc_id path_arc = newCol->numPathArcs;
   PathArcListNode* listNode = &newCol->pathArcs[path_arc];
   listNode->arc = arc;

   listNode->nextMember = newCol->reducedMembers[reducedMember].firstPathArc;
   newCol->reducedMembers[reducedMember].firstPathArc = path_arc;
   newCol->reducedMembers[reducedMember].numPathArcs += 1;

   listNode->nextOverall = newCol->firstOverallPathArc;
   newCol->firstOverallPathArc = path_arc;

   ++newCol->numPathArcs;
   assert(newCol->numPathArcs <= newCol->memPathArcs);

   assert(arc < newCol->memArcsInPath);
   newCol->arcInPath[arc] = TRUE;
   newCol->arcInPathReversed[arc] = reversed;
   assert(memberIsRepresentative(dec, newCol->reducedMembers[reducedMember].member));
   if( getMemberType(dec, newCol->reducedMembers[reducedMember].member) == SPQR_MEMBERTYPE_RIGID )
   {

      listNode->arcHead = findEffectiveArcHead(dec, arc);
      listNode->arcTail = findEffectiveArcTail(dec, arc);
      if( reversed )
      {
         SCIPswapInts(&listNode->arcHead, &listNode->arcTail);
      }
      assert(SPQRnodeIsValid(listNode->arcHead) && SPQRnodeIsValid(listNode->arcTail));
      assert(listNode->arcHead < newCol->memNodePathDegree && listNode->arcTail < newCol->memNodePathDegree);
      ++newCol->nodeInPathDegree[listNode->arcHead];
      ++newCol->nodeOutPathDegree[listNode->arcTail];
   } else
   {
      listNode->arcHead = SPQR_INVALID_NODE;
      listNode->arcTail = SPQR_INVALID_NODE;
   }
   listNode->reversed = reversed;
}

static SCIP_RETCODE createPathArcs(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol
)
{
   int maxNumPathArcs = newCol->numDecompositionRowArcs + getNumMembers(dec);
   if( newCol->memPathArcs < maxNumPathArcs )
   {
      int newMaxNumArcs = 2 * maxNumPathArcs;//safety factor to prevent very frequent reallocations
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newCol->pathArcs, newCol->memPathArcs, newMaxNumArcs));
      newCol->memPathArcs = newMaxNumArcs;
   }
   int maxPathArcIndex = largestArcID(dec);
   if( newCol->memArcsInPath < maxPathArcIndex )
   {
      int newSize = 2 * maxPathArcIndex;//safety factor to prevent very frequent reallocations
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newCol->arcInPath, newCol->memArcsInPath, newSize));
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newCol->arcInPathReversed, newCol->memArcsInPath, newSize));

      for( int i = newCol->memArcsInPath; i < newSize; ++i )
      {
         newCol->arcInPath[i] = FALSE;
         newCol->arcInPathReversed[i] = FALSE;
      }
      newCol->memArcsInPath = newSize;
   }
   int maxNumNodes = largestNodeID(dec);
   if( newCol->memNodePathDegree < maxNumNodes )
   {
      int newSize = 2 * maxNumNodes;//safety factor to prevent very frequent reallocations
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newCol->nodeInPathDegree, newCol->memNodePathDegree, newSize));
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newCol->nodeOutPathDegree, newCol->memNodePathDegree, newSize));
      for( int i = newCol->memNodePathDegree; i < newSize; ++i )
      {
         newCol->nodeInPathDegree[i] = 0;
         newCol->nodeOutPathDegree[i] = 0;
      }
      newCol->memNodePathDegree = newSize;
   }
   for( int i = 0; i < newCol->numDecompositionRowArcs; ++i )
   {
      spqr_arc arc = newCol->decompositionRowArcs[i];
      spqr_member member = findArcMember(dec, arc);
      reduced_member_id reducedMember = newCol->memberInformation[member].reducedMember;
      createPathArc(dec, newCol, arc, reducedMember, newCol->decompositionArcReversed[i]);
   }

   return SCIP_OKAY;
}


/**
 * Saves the information of the current row and partitions it based on whether or not the given columns are
 * already part of the decomposition.
 */
static SCIP_RETCODE
newColUpdateColInformation(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   spqr_col column,
   const spqr_row* nonzeroRows,
   const double* nonzeroValues,
   int numNonzeros
)
{
   newCol->newColIndex = column;

   newCol->numDecompositionRowArcs = 0;
   newCol->numNewRowArcs = 0;

   for( int i = 0; i < numNonzeros; ++i )
   {
      spqr_arc rowArc = getDecompositionRowArc(dec, nonzeroRows[i]);
      SCIP_Bool reversed = nonzeroValues[i] < 0.0;
      if( SPQRarcIsValid(rowArc))
      {//If the arc is the current decomposition: save it in the array
         if( newCol->numDecompositionRowArcs == newCol->memDecompositionRowArcs )
         {
            int newNumArcs = newCol->memDecompositionRowArcs == 0 ? 8 : 2 * newCol->memDecompositionRowArcs;
            SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newCol->decompositionRowArcs,
                                                  newCol->memDecompositionRowArcs, newNumArcs));
            SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newCol->decompositionArcReversed,
                                                   newCol->memDecompositionRowArcs, newNumArcs));
            newCol->memDecompositionRowArcs = newNumArcs;
         }
         newCol->decompositionRowArcs[newCol->numDecompositionRowArcs] = rowArc;
         newCol->decompositionArcReversed[newCol->numDecompositionRowArcs] = reversed;
         ++newCol->numDecompositionRowArcs;
      } else
      {
         //Not in the decomposition: add it to the set of arcs which are newly added with this row.
         if( newCol->numNewRowArcs == newCol->memNewRowArcs )
         {
            int newNumArcs = newCol->memNewRowArcs == 0 ? 8 : 2 * newCol->memNewRowArcs;
            SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newCol->newRowArcs,
                                                  newCol->memNewRowArcs, newNumArcs));
            SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newCol->newRowArcReversed,
                                                  newCol->memNewRowArcs, newNumArcs));
            newCol->memNewRowArcs = newNumArcs;
         }
         newCol->newRowArcs[newCol->numNewRowArcs] = nonzeroRows[i];
         newCol->newRowArcReversed[newCol->numNewRowArcs] = reversed;
         newCol->numNewRowArcs++;
      }
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE computeLeafMembers(
   const SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol
)
{
   if( newCol->numReducedMembers > newCol->memLeafMembers )
   {
      int newSize = maxValue(newCol->numReducedMembers, 2 * newCol->memLeafMembers);
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newCol->leafMembers, newCol->memLeafMembers, newSize));
      newCol->memLeafMembers = newSize;
   }
   newCol->numLeafMembers = 0;

   for( reduced_member_id reducedMember = 0; reducedMember < newCol->numReducedMembers; ++reducedMember )
   {
      if( newCol->reducedMembers[reducedMember].numChildren == 0 )
      {
         newCol->leafMembers[newCol->numLeafMembers] = reducedMember;
         ++newCol->numLeafMembers;
      }
   }
   return SCIP_OKAY;
}

static void determineRigidPath(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   SPQRColReducedMember* redMem
)
{
   assert(dec);
   assert(newCol);
   assert(redMem);

   SCIP_Bool isValidPath = TRUE;
   redMem->rigidPathStart = SPQR_INVALID_NODE;
   redMem->rigidPathEnd = SPQR_INVALID_NODE;
   for( path_arc_id pathArc = redMem->firstPathArc; pathArcIsValid(
      pathArc); pathArc = newCol->pathArcs[pathArc].nextMember )
   {
      spqr_node head = newCol->pathArcs[pathArc].arcHead;
      spqr_node tail = newCol->pathArcs[pathArc].arcTail;
      assert(nodeIsRepresentative(dec, head) && nodeIsRepresentative(dec, tail));

      if( newCol->nodeInPathDegree[head] > 1 || newCol->nodeOutPathDegree[tail] > 1 )
      {
         //not network -> stop
         isValidPath = FALSE;
         break;
      }
      if( newCol->nodeOutPathDegree[head] == 0 )
      {
         //found end node
         //If this is the second, stop
         if( SPQRnodeIsValid(redMem->rigidPathEnd))
         {
            isValidPath = FALSE;
            break;
         }
         redMem->rigidPathEnd = head;
      }
      if( newCol->nodeInPathDegree[tail] == 0 )
      {
         //Found start node.
         //If this is the second, stop.
         if( SPQRnodeIsValid(redMem->rigidPathStart))
         {
            isValidPath = FALSE;
            break;
         }
         redMem->rigidPathStart = tail;
      }
   }
   //assert that both a start and end node have been found
   assert(!isValidPath || ( SPQRnodeIsValid(redMem->rigidPathStart) && SPQRnodeIsValid(redMem->rigidPathEnd)));
   if( !isValidPath )
   {
      redMem->rigidPathStart = SPQR_INVALID_NODE;
      redMem->rigidPathEnd = SPQR_INVALID_NODE;
      redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
      newCol->remainsNetwork = FALSE;
   }
}

static void determineSingleRigidType(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   reduced_member_id reducedMember
   )
{
   assert(dec);
   assert(newCol);
   assert(reducedMemberIsValid(reducedMember));
   SPQRColReducedMember* redMem = &newCol->reducedMembers[reducedMember];
   assert(pathArcIsValid(redMem->firstPathArc));
   determineRigidPath(dec, newCol, redMem);
   if( redMem->type != REDUCEDMEMBER_TYPE_NOT_NETWORK )
   {
      redMem->type = REDUCEDMEMBER_TYPE_MERGED;
   }
}

//TODO: type seems somewhat duplicate
static void determineSingleComponentType(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   reduced_member_id reducedMember
)
{
   assert(dec);
   assert(newCol);

   int numNonPropagatedAdjacent =
      newCol->reducedMembers[reducedMember].numChildren - newCol->reducedMembers[reducedMember].numPropagatedChildren;
   if( reducedMemberIsValid(newCol->reducedMembers[reducedMember].parent) &&
       newCol->reducedMembers[newCol->reducedMembers[reducedMember].parent].type != REDUCEDMEMBER_TYPE_CYCLE )
   {
      ++numNonPropagatedAdjacent;
   }

   if( numNonPropagatedAdjacent > 2 )
   {
      newCol->reducedMembers[reducedMember].type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
      newCol->remainsNetwork = FALSE;
      return;
   }


   spqr_member member = findMember(dec, newCol->reducedMembers[reducedMember].member);
   SPQRMemberType type = getMemberType(dec, member);
   switch( type )
   {
      case SPQR_MEMBERTYPE_RIGID:
      {
         determineSingleRigidType(dec, newCol, reducedMember);
         break;
      }
      case SPQR_MEMBERTYPE_PARALLEL:
      {
         SPQRColReducedMember* redMem = &newCol->reducedMembers[reducedMember];
         assert(pathArcIsValid(redMem->firstPathArc));
         SCIP_Bool pathForwards = newCol->pathArcs[redMem->firstPathArc].reversed ==
                                  arcIsReversedNonRigid(dec, newCol->pathArcs[redMem->firstPathArc].arc);
         redMem->pathBackwards = !pathForwards;
         redMem->type = REDUCEDMEMBER_TYPE_CYCLE;
         break;
      }
      case SPQR_MEMBERTYPE_SERIES:
      case SPQR_MEMBERTYPE_LOOP:
      {
         SPQRColReducedMember* redMem = &newCol->reducedMembers[reducedMember];
         int countedPathArcs = 0;
         SCIP_Bool good = TRUE;
         SCIP_Bool passesForwards = TRUE;
         for( path_arc_id pathArc = redMem->firstPathArc; pathArcIsValid(pathArc);
              pathArc = newCol->pathArcs[pathArc].nextMember )
         {
            if( countedPathArcs == 0 )
            {
               passesForwards =
                  newCol->pathArcs[pathArc].reversed != arcIsReversedNonRigid(dec, newCol->pathArcs[pathArc].arc);
            } else if(
               ( newCol->pathArcs[pathArc].reversed != arcIsReversedNonRigid(dec, newCol->pathArcs[pathArc].arc)) !=
               passesForwards )
            {
               good = FALSE;
               break;
            }
            ++countedPathArcs;
         }
         if( !good )
         {
            redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
            newCol->remainsNetwork = FALSE;
            break;
         }

         redMem->pathBackwards = !passesForwards;
         if( countedPathArcs == getNumMemberArcs(dec, findMember(dec, redMem->member)) - 1 )
         {
            //Type -> Cycle;
            //Propagate arc
            redMem->type = REDUCEDMEMBER_TYPE_CYCLE;
         } else
         {
            //Type -> single_end
            redMem->type = REDUCEDMEMBER_TYPE_MERGED;
         }
         break;
      }
      case SPQR_MEMBERTYPE_UNASSIGNED:
      {
         assert(FALSE);
         break;
      }
   }
}


static void determinePathSeriesType(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   reduced_member_id reducedMember,
   spqr_member member,
   MemberPathType previousType,
   spqr_arc source,
   spqr_arc target
)
{
   assert(dec);
   assert(newCol);
   assert(reducedMemberIsValid(reducedMember));
   assert(SPQRmemberIsValid(member) && memberIsRepresentative(dec, member));
   assert(getMemberType(dec, member) == SPQR_MEMBERTYPE_SERIES);

   SPQRColReducedMember* redMem = &newCol->reducedMembers[reducedMember];
   int countedPathArcs = 0;

   SCIP_Bool good = TRUE;
   SCIP_Bool passesForwards = TRUE;
   for( path_arc_id pathArc = redMem->firstPathArc; pathArcIsValid(pathArc);
        pathArc = newCol->pathArcs[pathArc].nextMember )
   {
      if( countedPathArcs == 0 )
      {
         passesForwards =
            newCol->pathArcs[pathArc].reversed != arcIsReversedNonRigid(dec, newCol->pathArcs[pathArc].arc);
      } else if(( newCol->pathArcs[pathArc].reversed != arcIsReversedNonRigid(dec, newCol->pathArcs[pathArc].arc)) !=
                passesForwards )
      {
         good = FALSE;
         break;
      }
      ++countedPathArcs;
   }
   //If the internal directions of the arcs do not agree, we have no way to realize
   if( !good )
   {
      redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
      newCol->remainsNetwork = FALSE;
      return;
   }
   //If we are in the first skeleton processed, ignore the previous member type
   if( !SPQRarcIsValid(source))
   {
      assert(countedPathArcs > 0);
      assert(SPQRarcIsValid(target));
      SCIP_Bool firstReversed = arcIsReversedNonRigid(dec, newCol->pathArcs[redMem->firstPathArc].arc);
      SCIP_Bool targetReversed = arcIsReversedNonRigid(dec, target);
      SCIP_Bool reversePath = newCol->pathArcs[redMem->firstPathArc].reversed;

      if(( firstReversed == targetReversed ) == reversePath )
      {
         redMem->pathType = INTO_HEAD;
      } else
      {
         redMem->pathType = OUT_HEAD;
      }
      redMem->pathBackwards = !passesForwards;
      return;
   }
   SCIP_Bool sourceReversed = arcIsReversedNonRigid(dec, source);
   if( countedPathArcs > 0 )
   {
      SCIP_Bool isIntoHeadOrOutTail = isInto(previousType) == isHead(previousType);
      SCIP_Bool isGood = isIntoHeadOrOutTail == ( sourceReversed == passesForwards );
      if( !isGood )
      {
         redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
         newCol->remainsNetwork = FALSE;
         return;
      }
   }
   redMem->pathBackwards = !passesForwards;
   if( SPQRarcIsValid(target))
   {
      SCIP_Bool targetReversed = arcIsReversedNonRigid(dec, target);
      SCIP_Bool inSameDirection = sourceReversed == targetReversed;

      MemberPathType currentType;
      switch( previousType )
      {
         case INTO_HEAD:
         {
            currentType = inSameDirection ? INTO_TAIL : INTO_HEAD;
            break;
         }
         case INTO_TAIL:
         {
            currentType = inSameDirection ? INTO_HEAD : INTO_TAIL;
            break;
         }
         case OUT_HEAD:
         {
            currentType = inSameDirection ? OUT_TAIL : OUT_HEAD;
            break;
         }
         case OUT_TAIL:
         {
            currentType = inSameDirection ? OUT_HEAD : OUT_TAIL;
            break;
         }
         default:
         {
            assert(FALSE);
         }
      }
      redMem->pathType = currentType;
      return;
   }
   //If we are in the last skeleton, we only have a source, so nothing further to do
   assert(countedPathArcs > 0);
   //Strictly speaking below are no-ops within the algorithm, but help with debugging
   if( isInto(previousType))
   {
      redMem->pathType = INTO_HEAD;
   } else
   {
      redMem->pathType = OUT_HEAD;
   }
}

static void determinePathParallelType(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   reduced_member_id reducedMember,
   spqr_member member,
   MemberPathType previousType,
   spqr_arc source,
   spqr_arc target
)
{
   assert(dec);
   assert(newCol);
   assert(reducedMemberIsValid(reducedMember));
   assert(SPQRmemberIsValid(member) && memberIsRepresentative(dec, member));
   assert(getMemberType(dec, member) == SPQR_MEMBERTYPE_PARALLEL);

   //Parallel members must always be of degree two; if they are a leaf, they must contain an edge, which could have been propagated
   assert(SPQRarcIsValid(source) && SPQRarcIsValid(target));
   SCIP_Bool sourceReversed = arcIsReversedNonRigid(dec, source);
   SCIP_Bool targetReversed = arcIsReversedNonRigid(dec, target);
   SCIP_Bool inSameDirection = sourceReversed == targetReversed;

   SPQRColReducedMember* redMem = &newCol->reducedMembers[reducedMember];

   path_arc_id pathArc = redMem->firstPathArc;

   SCIP_Bool arcPresent = FALSE;
   if( pathArcIsValid(pathArc))
   {
      SCIP_Bool pathArcReversed =
         newCol->pathArcs[pathArc].reversed != arcIsReversedNonRigid(dec, newCol->pathArcs[pathArc].arc);
      SCIP_Bool intoHeadOrOutTail = isInto(previousType) == isHead(previousType);
      SCIP_Bool good = intoHeadOrOutTail == ( pathArcReversed != sourceReversed );
      if( !good )
      {
         redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
         newCol->remainsNetwork = FALSE;
         return;
      }
      arcPresent = TRUE;
   }

   SCIP_Bool swapHeadTail = arcPresent != inSameDirection;
   switch( previousType )
   {
      case INTO_HEAD:
      {
         redMem->pathType = swapHeadTail ? INTO_HEAD : INTO_TAIL;
         break;
      }
      case INTO_TAIL:
      {
         redMem->pathType = swapHeadTail ? INTO_TAIL : INTO_HEAD;
         break;
      }
      case OUT_HEAD:
      {
         redMem->pathType = swapHeadTail ? OUT_HEAD : OUT_TAIL;
         break;
      }
      case OUT_TAIL:
      {
         redMem->pathType = swapHeadTail ? OUT_TAIL : OUT_HEAD;
         break;
      }
   }
}

static void determinePathRigidType(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   reduced_member_id reducedMember,
   MemberPathType previousType,
   spqr_arc source,
   spqr_arc target
)
{

   SPQRColReducedMember* redMem = &newCol->reducedMembers[reducedMember];
   if( pathArcIsInvalid(redMem->firstPathArc))
   {
      assert(SPQRarcIsValid(source));
      assert(SPQRarcIsValid(target));
      //In this case, we need to check if the source and target are adjacent in any node

      spqr_node sourceTail = findEffectiveArcTail(dec, source);
      spqr_node sourceHead = findEffectiveArcHead(dec, source);
      spqr_node targetTail = findEffectiveArcTail(dec, target);
      spqr_node targetHead = findEffectiveArcHead(dec, target);
      SCIP_Bool sourceHeadIsTargetHead = sourceHead == targetHead;
      SCIP_Bool sourceTailIsTargetHead = sourceTail == targetHead;
      SCIP_Bool sourceHeadIsTargetTail = sourceHead == targetTail;
      SCIP_Bool sourceTailIsTargetTail = sourceTail == targetTail;

      if( !( sourceHeadIsTargetHead || sourceTailIsTargetHead || sourceHeadIsTargetTail || sourceTailIsTargetTail ))
      {
         redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
         newCol->remainsNetwork = FALSE;
         return;
      }
      assert(
         ( sourceHeadIsTargetHead ? 1 : 0 ) + ( sourceHeadIsTargetTail ? 1 : 0 ) + ( sourceTailIsTargetHead ? 1 : 0 ) +
         ( sourceTailIsTargetTail ? 1 : 0 ) == 1);
      //assert simplicity; they can be adjacent in only exactly one node
      SCIP_Bool isSourceHead = sourceHeadIsTargetHead || sourceHeadIsTargetTail;
      SCIP_Bool isTargetTail = sourceHeadIsTargetTail || sourceTailIsTargetTail;
      if( isHead(previousType) == isSourceHead )
      {
         //no need to reflect
         redMem->reverseArcs = FALSE;
         if( isInto(previousType))
         {
            if( isTargetTail )
            {
               redMem->pathType = INTO_TAIL;
            } else
            {
               redMem->pathType = INTO_HEAD;
            }
         } else
         {
            if( isTargetTail )
            {
               redMem->pathType = OUT_TAIL;
            } else
            {
               redMem->pathType = OUT_HEAD;
            }
         }
      } else
      {
         redMem->reverseArcs = TRUE;
         //because of the reversal, all the heads/tails are switched below
         if( isInto(previousType))
         {
            if( isTargetTail )
            {
               redMem->pathType = INTO_HEAD;
            } else
            {
               redMem->pathType = INTO_TAIL;
            }
         } else
         {
            if( isTargetTail )
            {
               redMem->pathType = OUT_HEAD;
            } else
            {
               redMem->pathType = OUT_TAIL;
            }
         }
      }
      assert(isInto(redMem->pathType) == isInto(previousType));
      return;
   }
   determineRigidPath(dec, newCol, redMem);
   if( redMem->type == REDUCEDMEMBER_TYPE_NOT_NETWORK )
   {
      return;
   }
   if( SPQRarcIsInvalid(source))
   {
      assert(SPQRarcIsValid(target));
      spqr_node targetTail = findEffectiveArcTail(dec, target);
      spqr_node targetHead = findEffectiveArcHead(dec, target);
      redMem->reverseArcs = FALSE;
      if( redMem->rigidPathEnd == targetHead )
      {
         redMem->pathType = INTO_HEAD;
      } else if( redMem->rigidPathEnd == targetTail )
      {
         redMem->pathType = INTO_TAIL;
      } else if( redMem->rigidPathStart == targetHead )
      {
         redMem->pathType = OUT_HEAD;
      } else if( redMem->rigidPathStart == targetTail )
      {
         redMem->pathType = OUT_TAIL;
      } else
      {
         redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
         newCol->remainsNetwork = FALSE;
      }
      return;
   }
   assert(SPQRarcIsValid(source));
   spqr_node sourceTail = findEffectiveArcTail(dec, source);
   spqr_node sourceHead = findEffectiveArcHead(dec, source);

   SCIP_Bool startsAtHead = sourceHead == redMem->rigidPathStart;
   SCIP_Bool endsAtTail = sourceTail == redMem->rigidPathEnd;
   SCIP_Bool startsAtTail = sourceTail == redMem->rigidPathStart;
   SCIP_Bool endsAtHead = sourceHead == redMem->rigidPathEnd;

   SCIP_Bool isIntoHeadOrOutTail = isInto(previousType) == isHead(previousType);
   if( isIntoHeadOrOutTail )
   {//into head or outTail
      //Check if path starts at head or ends at tail
      if( !startsAtHead && !endsAtTail )
      {
         redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
         newCol->remainsNetwork = FALSE;
         return;
      }
      assert(startsAtHead || endsAtTail);//both can hold; they can form cycle but other components can not be reduced
      //        redMem->reverseArcs = isInto(previousType) != startsAtHead; //Reverse only if there is no path starting at head
   } else
   {//Into tail or outHead
      //Check if path starts at tail or ends at head
      if( !startsAtTail && !endsAtHead )
      {
         redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
         newCol->remainsNetwork = FALSE;
         return;
      }
      assert(startsAtTail || endsAtHead);// both can hold; they can form cycle but other components can not be reduced
      //        redMem->reverseArcs = isInto(previousType) != startsAtTail; //Reverse only if there is no path starting at tail
   }

   if( SPQRarcIsValid(target))
   {
      spqr_node targetTail = findEffectiveArcTail(dec, target);
      spqr_node targetHead = findEffectiveArcHead(dec, target);

      //Check if they are not parallel; (below logic relies on this fact)
      assert(!(( targetHead == sourceHead && targetTail == sourceTail ) ||
               ( targetHead == sourceTail && targetTail == sourceHead )));

      SCIP_Bool startsAtTargetHead = redMem->rigidPathStart == targetHead;
      SCIP_Bool startsAtTargetTail = redMem->rigidPathStart == targetTail;
      SCIP_Bool endsAtTargetHead = redMem->rigidPathEnd == targetHead;
      SCIP_Bool endsAtTargetTail = redMem->rigidPathEnd == targetTail;

      if( !( startsAtTargetHead || startsAtTargetTail || endsAtTargetHead || endsAtTargetTail ))
      {
         redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
         newCol->remainsNetwork = FALSE;
         return;
      }
      SCIP_Bool outReverse = FALSE;
      SCIP_Bool outHead = FALSE;

      if( isInto(previousType) == isHead(previousType))
      {
         if( startsAtHead && endsAtTail )
         {
            outReverse = ( startsAtTargetHead || startsAtTargetTail ) == isInto(previousType);
         } else if( startsAtHead )
         {
            outReverse = !isInto(previousType);
         } else
         {
            assert(endsAtTail);
            outReverse = isInto(previousType);
         }
      } else
      {
         if( startsAtTail && endsAtHead )
         {
            outReverse = ( startsAtTargetHead || startsAtTargetTail ) == isInto(previousType);
         } else if( startsAtTail )
         {
            outReverse = !isInto(previousType);
         } else
         {
            assert(endsAtHead);
            outReverse = isInto(previousType);
         }
      }

      //TODO: these can probably be simplified significantly, but this might pose risk of introducing incorrect assumptions
      SCIP_Bool isBad = FALSE;
      if( isInto(previousType) == isHead(previousType))
      {
         if( startsAtHead && endsAtTail )
         {
            outHead = ( startsAtTargetTail || endsAtTargetHead ) == isInto(previousType);
         } else if( startsAtHead )
         {
            if( endsAtTargetHead )
            {
               outHead = isInto(previousType);
            } else if( endsAtTargetTail )
            {
               outHead = !isInto(previousType);
            } else
            {
               isBad = TRUE;
            }
         } else
         {
            assert(endsAtTail);
            if( startsAtTargetTail )
            {
               outHead = isInto(previousType);
            } else if( startsAtTargetHead )
            {
               outHead = !isInto(previousType);
            } else
            {
               isBad = TRUE;
            }
         }
      } else
      {
         if( startsAtTail && endsAtHead )
         {
            outHead = ( startsAtTargetTail || endsAtTargetHead ) == isInto(previousType);
         } else if( startsAtTail )
         {
            if( endsAtTargetHead )
            {
               outHead = isInto(previousType);
            } else if( endsAtTargetTail )
            {
               outHead = !isInto(previousType);
            } else
            {
               isBad = TRUE;
            }
         } else
         {
            assert(endsAtHead);
            if( startsAtTargetTail )
            {
               outHead = isInto(previousType);
            } else if( startsAtTargetHead )
            {
               outHead = !isInto(previousType);
            } else
            {
               isBad = TRUE;
            }
         }
      }
      if( isBad )
      {
         redMem->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
         newCol->remainsNetwork = FALSE;
         return;
      }

      redMem->reverseArcs = outReverse;
      if( isInto(previousType))
      {
         redMem->pathType = outHead ? INTO_HEAD : INTO_TAIL;
      } else
      {
         redMem->pathType = outHead ? OUT_HEAD : OUT_TAIL;
      }
      return;
   }

   //TODO: is this simplifyable?
   if( isInto(previousType) == isHead(previousType))
   {
      redMem->reverseArcs = startsAtHead != isInto(previousType);
   } else
   {
      redMem->reverseArcs = startsAtTail != isInto(previousType);
   }
   //last member of the path. Since we already checked the source,
   //Below is technically no op, but helps with debugging
   if( isInto(previousType))
   {
      redMem->pathType = INTO_HEAD;
   } else
   {
      redMem->pathType = OUT_HEAD;
   }
}

static void determinePathMemberType(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   reduced_member_id reducedMember,
   spqr_member member,
   MemberPathType previousType,
   spqr_arc source,
   spqr_arc target
)
{
   newCol->reducedMembers[reducedMember].pathSourceArc = source;
   newCol->reducedMembers[reducedMember].pathTargetArc = target;
   //Check if the marked edges with the given signs
   //form a (reverse) directed path from one of the source's end nodes to one of the target's end nodes
   switch( getMemberType(dec, member))
   {
      case SPQR_MEMBERTYPE_RIGID:
      {
         determinePathRigidType(dec, newCol, reducedMember, previousType, source, target);
         break;
      }
      case SPQR_MEMBERTYPE_PARALLEL:
      {
         determinePathParallelType(dec, newCol, reducedMember, member, previousType, source, target);
         break;
      }
      case SPQR_MEMBERTYPE_SERIES:
      {
         determinePathSeriesType(dec, newCol, reducedMember, member, previousType, source, target);
         break;
      }
      case SPQR_MEMBERTYPE_LOOP:
      case SPQR_MEMBERTYPE_UNASSIGNED:
      {
         //In release
         newCol->remainsNetwork = FALSE;
         assert(FALSE);
         break;
      }
   }
}

static void determinePathTypes(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   SPQRColReducedComponent* component
)
{
   assert(dec);
   assert(newCol);
   assert(component);
   assert(component->numPathEndMembers == 2);

   assert(component->pathEndMembers[0] != component->root);

   //We check the path by going from end to end. We start at the leaf and process every component,
   //walking down until we hit the root.
   //Then, we walk up from the root and process each component in the same manner.
   reduced_member_id reducedStart = component->pathEndMembers[0];
   spqr_arc toPrevious = SPQR_INVALID_ARC;
   spqr_arc toNext = SPQR_INVALID_ARC;
   MemberPathType previousType = INTO_HEAD;//Arbitrary, is ignored in the first call

   spqr_member member = newCol->reducedMembers[reducedStart].member;
   reduced_member_id reducedMember = reducedStart;
   reduced_member_id previousReducedMember = reducedStart;
   while( reducedMember != component->root )
   {
      toNext = markerToParent(dec, member);
      determinePathMemberType(dec, newCol, reducedMember, member, previousType, toPrevious, toNext);
      if( !newCol->remainsNetwork )
      {
         return;
      }
      previousType = newCol->reducedMembers[reducedMember].pathType;
      toPrevious = markerOfParent(dec, member);
      member = findMemberParent(dec, newCol->reducedMembers[reducedMember].member);
      previousReducedMember = reducedMember;
      reducedMember = newCol->memberInformation[member].reducedMember;
      newCol->reducedMembers[previousReducedMember].nextPathMember = reducedMember;
      newCol->reducedMembers[previousReducedMember].nextPathMemberIsParent = TRUE;
   }

   while( reducedMember != component->pathEndMembers[1] )
   {
      //Search the (other) child node
      reduced_member_id child = INVALID_REDUCED_MEMBER;
      //a bit ugly linear search, but not a problem for time complexity
      for( children_idx i = newCol->reducedMembers[reducedMember].firstChild;
           i <
           newCol->reducedMembers[reducedMember].firstChild + newCol->reducedMembers[reducedMember].numChildren; ++i )
      {
         reduced_member_id childReduced = newCol->childrenStorage[i];
         if( newCol->reducedMembers[childReduced].type != REDUCEDMEMBER_TYPE_CYCLE &&
             childReduced != previousReducedMember )
         {
            child = childReduced;
            break;
         }
      }
      assert(reducedMemberIsValid(child));

      spqr_member childMember = newCol->reducedMembers[child].member;
      toNext = markerOfParent(dec, childMember);

      determinePathMemberType(dec, newCol, reducedMember, member, previousType, toPrevious, toNext);
      if( !newCol->remainsNetwork )
      {
         return;
      }
      previousType = newCol->reducedMembers[reducedMember].pathType;
      toPrevious = markerToParent(dec, childMember);
      member = childMember;
      previousReducedMember = reducedMember;
      reducedMember = child;
      newCol->reducedMembers[previousReducedMember].nextPathMember = reducedMember;
      newCol->reducedMembers[previousReducedMember].nextPathMemberIsParent = FALSE;
   }
   //The last iteration is not performed by the loops above.
   //We explicitly set the target arc to invalid in order to indicate that this is the last iteration.
   toNext = SPQR_INVALID_ARC;
   determinePathMemberType(dec, newCol, reducedMember, member, previousType, toPrevious, toNext);
   newCol->reducedMembers[reducedMember].nextPathMember = INVALID_REDUCED_MEMBER;
   //since we return anyways, no need to check newCol->remainsNetwork explicitly
}

static void checkRigidLeaf(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   reduced_member_id leaf,
   spqr_arc toParent,
   reduced_member_id parent,
   spqr_arc toChild
)
{
   SPQRColReducedMember* leafMember = &newCol->reducedMembers[leaf];
   determineRigidPath(dec, newCol, leafMember);
   if( leafMember->type == REDUCEDMEMBER_TYPE_NOT_NETWORK )
   {
      return;
   }
   spqr_node targetHead = findEffectiveArcHead(dec, toParent);
   spqr_node targetTail = findEffectiveArcTail(dec, toParent);
   SCIP_Bool matches = leafMember->rigidPathStart == targetTail && leafMember->rigidPathEnd == targetHead;
   SCIP_Bool opposite = leafMember->rigidPathStart == targetHead && leafMember->rigidPathEnd == targetTail;
   if( matches || opposite )
   {
      leafMember->type = REDUCEDMEMBER_TYPE_CYCLE;
      createPathArc(dec, newCol, toChild, parent, opposite);
      return;
   }
   leafMember->type = REDUCEDMEMBER_TYPE_MERGED;
}

static ReducedMemberType checkLeaf(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   reduced_member_id leaf,
   spqr_arc toParent,
   reduced_member_id parent,
   spqr_arc toChild
)
{
   assert(dec);
   assert(newCol);
   assert(SPQRarcIsValid(toParent));
   assert(SPQRarcIsValid(toChild));
   assert(!arcIsTree(dec, toParent));
   assert(reducedMemberIsValid(leaf));
   assert(reducedMemberIsValid(parent));

   switch( getMemberType(dec, newCol->reducedMembers[leaf].member))
   {
      case SPQR_MEMBERTYPE_RIGID:
      {
         checkRigidLeaf(dec, newCol, leaf, toParent, parent, toChild);
         break;
      }
      case SPQR_MEMBERTYPE_PARALLEL:
      {
         SPQRColReducedMember* reducedMember = &newCol->reducedMembers[leaf];
         assert(pathArcIsValid(reducedMember->firstPathArc));
         reducedMember->type = REDUCEDMEMBER_TYPE_CYCLE;

         SCIP_Bool pathArcReversed = newCol->pathArcs[reducedMember->firstPathArc].reversed;
         SCIP_Bool arcPathArcIsReverse = arcIsReversedNonRigid(dec, newCol->pathArcs[reducedMember->firstPathArc].arc);
         SCIP_Bool parentReversed = arcIsReversedNonRigid(dec, toParent);
         createPathArc(dec, newCol, toChild, parent, ( arcPathArcIsReverse == parentReversed ) == pathArcReversed);
         break;
      }
      case SPQR_MEMBERTYPE_SERIES:
      case SPQR_MEMBERTYPE_LOOP:
      {
         SPQRColReducedMember* reducedMember = &newCol->reducedMembers[leaf];
         int countedPathArcs = 0;
         SCIP_Bool good = TRUE;
         SCIP_Bool passesForwards = TRUE;
         for( path_arc_id pathArc = reducedMember->firstPathArc; pathArcIsValid(pathArc);
              pathArc = newCol->pathArcs[pathArc].nextMember )
         {
            if( countedPathArcs == 0 )
            {
               passesForwards =
                  newCol->pathArcs[pathArc].reversed != arcIsReversedNonRigid(dec, newCol->pathArcs[pathArc].arc);
            } else if(
               ( newCol->pathArcs[pathArc].reversed != arcIsReversedNonRigid(dec, newCol->pathArcs[pathArc].arc)) !=
               passesForwards )
            {
               good = FALSE;
               break;
            }
            ++countedPathArcs;
         }
         if( !good )
         {
            reducedMember->type = REDUCEDMEMBER_TYPE_NOT_NETWORK;
            newCol->remainsNetwork = FALSE;
            break;
         }

         reducedMember->pathBackwards = !passesForwards;
         if( countedPathArcs == getNumMemberArcs(dec, findMember(dec, reducedMember->member)) - 1 )
         {
            //Type -> Cycle;
            //Propagate arc
            reducedMember->type = REDUCEDMEMBER_TYPE_CYCLE;

            SCIP_Bool firstArcReversed = arcIsReversedNonRigid(dec, newCol->pathArcs[reducedMember->firstPathArc].arc);
            SCIP_Bool firstArcInPathReverse = newCol->pathArcs[reducedMember->firstPathArc].reversed;
            SCIP_Bool parentReversed = arcIsReversedNonRigid(dec, toParent);
            createPathArc(dec, newCol, toChild, parent,
                          ( parentReversed == firstArcReversed ) != firstArcInPathReverse);

         } else
         {
            //Type -> single_end
            reducedMember->type = REDUCEDMEMBER_TYPE_MERGED;
         }

         break;
      }
      case SPQR_MEMBERTYPE_UNASSIGNED:
      {
         assert(FALSE);
         break;
      }
   }
   return newCol->reducedMembers[leaf].type;
}

static void propagateCycles(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol
)
{
   assert(dec);
   assert(newCol);
   int leafArrayIndex = 0;

   while( leafArrayIndex != newCol->numLeafMembers )
   {
      reduced_member_id leaf = newCol->leafMembers[leafArrayIndex];
      //next is invalid if the member is not in the reduced decomposition.
      reduced_member_id next = newCol->reducedMembers[leaf].parent;
      spqr_arc parentMarker = markerToParent(dec, newCol->reducedMembers[leaf].member);
      if( reducedMemberIsValid(next) && !arcIsTree(dec, parentMarker))
      {
         assert(reducedMemberIsValid(next));
         assert(SPQRarcIsValid(parentMarker));
         ReducedMemberType type = checkLeaf(dec, newCol, leaf, parentMarker, next,
                                            markerOfParent(dec, newCol->reducedMembers[leaf].member));
         if( type == REDUCEDMEMBER_TYPE_CYCLE )
         {
            ++newCol->reducedMembers[next].numPropagatedChildren;
            if( newCol->reducedMembers[next].numPropagatedChildren == newCol->reducedMembers[next].numChildren )
            {
               newCol->leafMembers[leafArrayIndex] = next;
            } else
            {
               ++leafArrayIndex;
            }
         } else if( type == REDUCEDMEMBER_TYPE_NOT_NETWORK )
         {
            return;
         } else
         {
            assert(type == REDUCEDMEMBER_TYPE_MERGED);
            ++leafArrayIndex;
            int component = newCol->reducedMembers[leaf].componentIndex;
            if( newCol->reducedComponents[component].numPathEndMembers >= 2 )
            {
               newCol->remainsNetwork = FALSE;
               return;
            }
            assert(newCol->reducedComponents[component].numPathEndMembers < 2);
            newCol->reducedComponents[component].pathEndMembers[newCol->reducedComponents[component].numPathEndMembers] = leaf;
            ++newCol->reducedComponents[component].numPathEndMembers;
         }
      } else
      {
         ++leafArrayIndex;
         int component = newCol->reducedMembers[leaf].componentIndex;
         if( newCol->reducedComponents[component].numPathEndMembers >= 2 )
         {
            newCol->remainsNetwork = FALSE;
            return;
         }
         assert(newCol->reducedComponents[component].numPathEndMembers < 2);
         newCol->reducedComponents[component].pathEndMembers[newCol->reducedComponents[component].numPathEndMembers] = leaf;
         ++newCol->reducedComponents[component].numPathEndMembers;
      }
   }

   for( int j = 0; j < newCol->numReducedComponents; ++j )
   {
      //The reduced root might be a leaf as well: we propagate it last
      reduced_member_id root = newCol->reducedComponents[j].root;

      while( TRUE ) /*lint !e716*/
      {
         if( newCol->reducedMembers[root].numPropagatedChildren != newCol->reducedMembers[root].numChildren - 1 )
         {
            break;
         }
         //TODO: bit ugly, have to do a linear search for the child
         reduced_member_id child = INVALID_REDUCED_MEMBER;
         spqr_arc markerToChild = SPQR_INVALID_ARC;
         for( children_idx i = newCol->reducedMembers[root].firstChild;
              i < newCol->reducedMembers[root].firstChild + newCol->reducedMembers[root].numChildren; ++i )
         {
            reduced_member_id childReduced = newCol->childrenStorage[i];
            if( newCol->reducedMembers[childReduced].type != REDUCEDMEMBER_TYPE_CYCLE )
            {
               child = childReduced;
               markerToChild = markerOfParent(dec, newCol->reducedMembers[child].member);
               break;
            }
         }
         assert(SPQRmemberIsValid(newCol->reducedMembers[child].member));
         assert(SPQRarcIsValid(markerToChild));
         if( !arcIsTree(dec, markerToChild))
         {
            ReducedMemberType type = checkLeaf(dec, newCol, root, markerToChild, child,
                                               markerToParent(dec, newCol->reducedMembers[child].member));
            if( type == REDUCEDMEMBER_TYPE_CYCLE )
            {
               root = child;
               continue;
            } else if( type == REDUCEDMEMBER_TYPE_NOT_NETWORK )
            {
               return;
            }
         }
         //If the root has exactly one neighbour and is not contained, it is also considered a path end member
         int component = newCol->reducedMembers[root].componentIndex;
         SCIP_Bool rootPresent = FALSE;
         for( int i = 0; i < newCol->reducedComponents[component].numPathEndMembers; ++i )
         {
            rootPresent = rootPresent || ( newCol->reducedComponents[component].pathEndMembers[i] == root );
         }
         if( !rootPresent )
         {
            if( newCol->reducedComponents[component].numPathEndMembers >= 2 )
            {
               newCol->remainsNetwork = FALSE;
               return;
            }
            newCol->reducedComponents[component].pathEndMembers[newCol->reducedComponents[component].numPathEndMembers] = root;
            ++newCol->reducedComponents[component].numPathEndMembers;
         }
         break;
      }

      newCol->reducedComponents[j].root = root;
      newCol->reducedMembers[root].parent = INVALID_REDUCED_MEMBER;
   }
}

static void determineComponentTypes(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   SPQRColReducedComponent* component
)
{
   assert(dec);
   assert(newCol);
   assert(component);

   if( component->numPathEndMembers == 1 )
   {
      assert(component->root == component->pathEndMembers[0]);
      determineSingleComponentType(dec, newCol, component->root);
   } else
   {
      assert(component->numPathEndMembers == 2);
      determinePathTypes(dec, newCol, component);
   }
}

static SCIP_RETCODE SCIPnetcoladdCheck(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* coladd,
   int column,
   const int* nonzrows,
   const double* nonzvals,
   int nnonzs
)
{
   assert(dec);
   assert(coladd);
   assert(nnonzs == 0 || ( nonzrows && nonzvals ));

   coladd->remainsNetwork = TRUE;
   cleanupPreviousIteration(dec, coladd);
   //assert that previous iteration was cleaned up

   //Store call data
   SCIP_CALL(newColUpdateColInformation(dec, coladd, column, nonzrows, nonzvals, nnonzs));

   //compute reduced decomposition
   SCIP_CALL(constructReducedDecomposition(dec, coladd));
   //initialize path arcs in reduced decomposition
   SCIP_CALL(createPathArcs(dec, coladd));
   SCIP_CALL(computeLeafMembers(dec, coladd));
   propagateCycles(dec, coladd);
   //determine types
   if( coladd->remainsNetwork )
   {
      for( int i = 0; i < coladd->numReducedComponents; ++i )
      {
         determineComponentTypes(dec, coladd, &coladd->reducedComponents[i]);
      }
   }
   //clean up memberInformation
   cleanUpMemberInformation(coladd);

   return SCIP_OKAY;
}

///Contains the data which tells us where to store the new column after the graph has been modified
///In case member is a parallel or series node, the respective new column and rows are placed in parallel (or series) with it
///Otherwise, the rigid member has a free spot between firstNode and secondNode
typedef struct
{
   spqr_member member;
   spqr_node head;
   spqr_node tail;
   spqr_arc representative;
   SCIP_Bool reversed;
} NewColInformation;

static NewColInformation emptyNewColInformation(void)
{
   NewColInformation information;
   information.member = SPQR_INVALID_MEMBER;
   information.head = SPQR_INVALID_NODE;
   information.tail = SPQR_INVALID_NODE;
   information.reversed = FALSE;
   information.representative = SPQR_INVALID_ARC;
   return information;
}

static void setTerminalHead(
   NewColInformation* info,
   spqr_node node
)
{
   assert(SPQRnodeIsValid(node));
   assert(SPQRnodeIsInvalid(info->head));
   assert(info);
   info->head = node;
}

static void setTerminalTail(
   NewColInformation* info,
   spqr_node node
)
{
   assert(SPQRnodeIsValid(node));
   assert(info);
   assert(SPQRnodeIsInvalid(info->tail));
   info->tail = node;
}

static void setTerminalReversed(
   NewColInformation* info,
   SCIP_Bool reversed
)
{
   assert(info);
   info->reversed = reversed;
}

static void setTerminalMember(
   NewColInformation* info,
   spqr_member member
)
{
   assert(info);
   info->member = member;
}

static void setTerminalRepresentative(
   NewColInformation* info,
   spqr_arc representative
)
{
   assert(info);
   info->representative = representative;
}

static SCIP_RETCODE splitParallel(
   SCIP_NETMATDECDATA* dec,
   spqr_member parallel,
   spqr_arc arc1,
   spqr_arc arc2,
   spqr_member* childParallel
)
{
   assert(dec);
   assert(SPQRarcIsValid(arc1));
   assert(SPQRarcIsValid(arc2));
   assert(SPQRmemberIsValid(parallel));

   SCIP_Bool childContainsTree = arcIsTree(dec, arc1) || arcIsTree(dec, arc2);
   spqr_arc toParent = markerToParent(dec, parallel);
   SCIP_Bool parentMoved = toParent == arc1 || toParent == arc2;
   SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, childParallel));

   moveArcToNewMember(dec, arc1, parallel, *childParallel);
   moveArcToNewMember(dec, arc2, parallel, *childParallel);

   if( parentMoved )
   {
      SCIP_CALL(createMarkerPair(dec, *childParallel, parallel, !childContainsTree, FALSE, FALSE));
   } else
   {
      SCIP_CALL(createMarkerPair(dec, parallel, *childParallel, childContainsTree, FALSE, FALSE));
   }
   return SCIP_OKAY;
}

static SCIP_RETCODE splitSeries(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   SPQRColReducedMember* reducedMember,
   spqr_member member,
   spqr_member* loopMember,
   NewColInformation* newColInfo
)
{
   assert(dec);
   assert(reducedMember);
   assert(SPQRmemberIsValid(member));
   assert(memberIsRepresentative(dec, member));

   assert(reducedMember->numPathArcs != 0);
   SCIP_Bool createPathSeries = reducedMember->numPathArcs > 1;
   SCIP_Bool convertOriginal = reducedMember->numPathArcs == getNumMemberArcs(dec, member) - 1;

   //TODO: for now very elaborate to first get logic correct. Can probably merge some branches later...
   if( !createPathSeries && !convertOriginal )
   {
      spqr_member parallel;
      SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &parallel));
      path_arc_id pathArcId = reducedMember->firstPathArc;
      assert(pathArcIsValid(pathArcId));
      spqr_arc marked = newCol->pathArcs[pathArcId].arc;
      assert(marked != markerToParent(dec, member));//TODO: handle this case later
      moveArcToNewMember(dec, marked, member, parallel);
      SCIP_Bool reversed = arcIsReversedNonRigid(dec, marked);
      SCIP_CALL(createMarkerPair(dec, member, parallel, TRUE, reversed, reversed));
      *loopMember = parallel;

      if( reversed == reducedMember->pathBackwards )
      {
         setTerminalReversed(newColInfo, !reversed);
      } else
      {
         setTerminalReversed(newColInfo, reversed);
      }

      setTerminalMember(newColInfo, *loopMember);
      return SCIP_OKAY;
   }
   if( !createPathSeries && convertOriginal )
   {
      //only one path arc; we are in a loop; no need to change anything
      assert(getNumMemberArcs(dec, member) == 2);
      assert(reducedMember->numPathArcs == 1);
      *loopMember = member;
      changeLoopToParallel(dec, member);

      path_arc_id pathArcId = reducedMember->firstPathArc;
      spqr_arc marked = newCol->pathArcs[pathArcId].arc;
      //The 'reversed' field has different meaning for parallels, so we need to change orientation when converting to a parallel
      arcFlipReversed(dec, marked);

      setTerminalReversed(newColInfo, reducedMember->pathBackwards);
      setTerminalMember(newColInfo, *loopMember);
      return SCIP_OKAY;
   }
   if( createPathSeries && !convertOriginal )
   {
      spqr_member pathMember;
      SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &pathMember));

      path_arc_id pathArcId = reducedMember->firstPathArc;
      SCIP_Bool parentMoved = FALSE;
      while( pathArcIsValid(pathArcId))
      {
         spqr_arc pathArc = newCol->pathArcs[pathArcId].arc;
         pathArcId = newCol->pathArcs[pathArcId].nextMember;
         if( pathArc == markerToParent(dec, member))
         {
            parentMoved = TRUE;
         }
         moveArcToNewMember(dec, pathArc, member, pathMember);
      }

      SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, loopMember));

      if( !parentMoved )
      {
         SCIP_CALL(createMarkerPair(dec, member, *loopMember, TRUE, FALSE, FALSE));
         SCIP_CALL(createMarkerPair(dec, *loopMember, pathMember, TRUE, FALSE, TRUE));
      } else
      {
         SCIP_CALL(createMarkerPair(dec, pathMember, *loopMember, FALSE, FALSE, TRUE));
         SCIP_CALL(createMarkerPair(dec, *loopMember, member, FALSE, FALSE, FALSE));
      }

      setTerminalReversed(newColInfo, !reducedMember->pathBackwards);
      setTerminalMember(newColInfo, *loopMember);
      return SCIP_OKAY;
   }
   assert(createPathSeries && convertOriginal);
   //There's one exception in this case
   //if the single unmarked (column) marker is a parent or child marker to a parallel member, we add the edge there
   {
      spqr_member adjacentMember = SPQR_INVALID_MEMBER;
      spqr_arc adjacentMarker = SPQR_INVALID_ARC;
      spqr_arc memberMarker = SPQR_INVALID_ARC;
      spqr_arc firstArc = getFirstMemberArc(dec, reducedMember->member);
      spqr_arc arc = firstArc;
      do
      {
         if( !newCol->arcInPath[arc] )
         {
            if( arc == markerToParent(dec, reducedMember->member))
            {
               adjacentMember = findMemberParent(dec, reducedMember->member);
               adjacentMarker = markerOfParent(dec, reducedMember->member);
               memberMarker = arc;
            } else if( arcIsMarker(dec, arc))
            {
               adjacentMember = findArcChildMember(dec, arc);
               adjacentMarker = markerToParent(dec, adjacentMember);
               memberMarker = arc;
            }

            break;//There is only a singular such arc
         }
         arc = getNextMemberArc(dec, arc);
      } while( arc != firstArc );

      if( SPQRmemberIsValid(adjacentMember))
      {
         SPQRMemberType adjacentType = getMemberType(dec, adjacentMember);
         if( adjacentType == SPQR_MEMBERTYPE_PARALLEL )
         {
            //Figure out if the markers are the same or opposite orientations
            //If they are the same, we can proceed as normal, otherwise, we need to flip the placed edge
            SCIP_Bool markersHaveSameOrientation =
               arcIsReversedNonRigid(dec, adjacentMarker) == arcIsReversedNonRigid(dec, memberMarker);
            setTerminalReversed(newColInfo, reducedMember->pathBackwards == markersHaveSameOrientation);
            setTerminalMember(newColInfo, adjacentMember);
            return SCIP_OKAY;
         }
      }
   }

   spqr_member pathMember;
   SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &pathMember));

   path_arc_id pathArcId = reducedMember->firstPathArc;
   SCIP_Bool parentMoved = FALSE;
   while( pathArcIsValid(pathArcId))
   {
      spqr_arc pathArc = newCol->pathArcs[pathArcId].arc;
      pathArcId = newCol->pathArcs[pathArcId].nextMember;
      if( pathArc == markerToParent(dec, member))
      {
         parentMoved = TRUE;
      }
      moveArcToNewMember(dec, pathArc, member, pathMember);
   }
   if( parentMoved )
   {
      SCIP_CALL(createMarkerPair(dec, pathMember, member, FALSE, FALSE, FALSE));
   } else
   {
      SCIP_CALL(createMarkerPair(dec, member, pathMember, TRUE, FALSE, FALSE));
   }

   changeLoopToParallel(dec, member);

   *loopMember = member;
   setTerminalReversed(newColInfo, reducedMember->pathBackwards);
   setTerminalMember(newColInfo, *loopMember);
   return SCIP_OKAY;
}


static SCIP_RETCODE splitSeriesMerging(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   SPQRColReducedMember* reducedMember,
   spqr_member member,
   spqr_arc* pathRepresentative,
   spqr_arc* nonPathRepresentative,
   spqr_arc exceptionArc1,
   spqr_arc exceptionArc2
)
{
   assert(dec);
   assert(reducedMember);
   assert(SPQRmemberIsValid(member));
   assert(memberIsRepresentative(dec, member));

   int numExceptionArcs = ( exceptionArc1 == SPQR_INVALID_ARC ? 0 : 1 ) + ( exceptionArc2 == SPQR_INVALID_ARC ? 0 : 1 );
   int numNonPathArcs = getNumMemberArcs(dec, member) - reducedMember->numPathArcs - numExceptionArcs;
   SCIP_Bool createPathSeries = reducedMember->numPathArcs > 1;
   //If this holds, there are 2 or more non-parent marker non-path arcs
   SCIP_Bool createNonPathSeries = numNonPathArcs > 1;
   assert(exceptionArc1 == SPQR_INVALID_ARC || !newCol->arcInPath[exceptionArc1]);
   assert(exceptionArc2 == SPQR_INVALID_ARC || !newCol->arcInPath[exceptionArc2]);

   if( createPathSeries )
   {
      spqr_member pathMember;
      SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &pathMember));

      path_arc_id pathArcId = reducedMember->firstPathArc;
      SCIP_Bool parentMoved = FALSE;
      while( pathArcIsValid(pathArcId))
      {
         spqr_arc pathArc = newCol->pathArcs[pathArcId].arc;
         pathArcId = newCol->pathArcs[pathArcId].nextMember;
         assert(pathArc != exceptionArc1 && pathArc != exceptionArc2);
         parentMoved = parentMoved || markerToParent(dec, member) == pathArc;
         moveArcToNewMember(dec, pathArc, member, pathMember);
      }
      assert(getNumMemberArcs(dec, pathMember) >= 2);

      spqr_arc ignored;
      SCIP_Bool inOldReversed = TRUE;
      SCIP_Bool inNewReversed = FALSE;
      if( parentMoved )
      {
         SCIP_CALL(createMarkerPairWithReferences(dec, pathMember, member, FALSE, inNewReversed, inOldReversed,
                                                  &ignored, pathRepresentative));
      } else
      {
         SCIP_CALL(createMarkerPairWithReferences(dec, member, pathMember, TRUE, inOldReversed, inNewReversed,
                                                  pathRepresentative, &ignored));
      }
   } else
   {
      if( pathArcIsValid(reducedMember->firstPathArc))
      {
         *pathRepresentative = newCol->pathArcs[reducedMember->firstPathArc].arc;
      } else
      {
         *pathRepresentative = SPQR_INVALID_ARC;
      }
   }

   if( createNonPathSeries )
   {
      spqr_member nonPathMember;
      SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &nonPathMember));

      spqr_arc arc = getFirstMemberArc(dec, member);
      SCIP_Bool parentMoved = FALSE;
      SCIP_Bool canStop = FALSE;//hack when the first arc is moved in the below loop to prevent that we immediately terminate
      do
      {
         spqr_arc nextArc = getNextMemberArc(dec, arc);
         if( arc != *pathRepresentative && arc != exceptionArc1 && arc != exceptionArc2 )
         {
            parentMoved = parentMoved || markerToParent(dec, member) == arc;
            moveArcToNewMember(dec, arc, member, nonPathMember);
         } else
         {
            canStop = TRUE;
         }
         arc = nextArc;
         if( canStop && arc == getFirstMemberArc(dec, member))
         {
            break;
         }
      } while( TRUE ); /*lint !e506*/
      assert(getNumMemberArcs(dec, nonPathMember) >= 2);
      SCIP_Bool representativeIsTree = !arcIsTree(dec, exceptionArc1);
      if( SPQRarcIsValid(exceptionArc2))
      {
         representativeIsTree = representativeIsTree || !arcIsTree(dec, exceptionArc2);
      }
      spqr_arc ignored;
      SCIP_Bool inOldReversed = TRUE;
      SCIP_Bool inNewReversed = FALSE;
      if( parentMoved )
      {
         SCIP_CALL(createMarkerPairWithReferences(dec, nonPathMember, member, !representativeIsTree,
                                                  inNewReversed, inOldReversed, &ignored, nonPathRepresentative));
      } else
      {
         SCIP_CALL(createMarkerPairWithReferences(dec, member, nonPathMember, representativeIsTree,
                                                  inOldReversed, inNewReversed, nonPathRepresentative, &ignored));
      }
   } else
   {
      *nonPathRepresentative = SPQR_INVALID_ARC;
      if( numNonPathArcs != 0 )
      {
         spqr_arc firstArc = getFirstMemberArc(dec, member);
         spqr_arc arc = firstArc;
         do
         {
            if( arc != *pathRepresentative && arc != exceptionArc1 && arc != exceptionArc2 )
            {
               *nonPathRepresentative = arc;
               break;
            }
            arc = getNextMemberArc(dec, arc);
         } while( arc != firstArc );
         assert(*nonPathRepresentative != SPQR_INVALID_ARC);
      }
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE transformFirstPathMember(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   reduced_member_id reducedMember,
   NewColInformation* newColInfo,
   spqr_arc* representativeArc,
   spqr_member* mergedMember
)
{
   spqr_member member = newCol->reducedMembers[reducedMember].member;
   SPQRMemberType type = getMemberType(dec, member);
   if( type == SPQR_MEMBERTYPE_RIGID )
   {
      //The nodes are already created, we only need to assign the correct start/end node
      switch( newCol->reducedMembers[reducedMember].pathType )
      {
         case INTO_HEAD:
         case INTO_TAIL:
            setTerminalTail(newColInfo, newCol->reducedMembers[reducedMember].rigidPathStart);
            break;
         case OUT_HEAD:
         case OUT_TAIL:
            setTerminalHead(newColInfo, newCol->reducedMembers[reducedMember].rigidPathEnd);
            break;
      }
      *representativeArc = findArcSign(dec, newCol->reducedMembers[reducedMember].pathTargetArc).representative;
      *mergedMember = member;

      return SCIP_OKAY;
   }
   assert(type == SPQR_MEMBERTYPE_SERIES);
   //Split off sets of multiple path non-path edges so that the series has exactly 3 edges

   spqr_arc target = newCol->reducedMembers[reducedMember].pathTargetArc;
   SPQRColReducedMember* redMem = &newCol->reducedMembers[reducedMember];
   spqr_arc pathRepresentative = SPQR_INVALID_ARC;
   spqr_arc nonPathRepresentative = SPQR_INVALID_ARC;
   SCIP_CALL(splitSeriesMerging(dec, newCol, redMem, member, &pathRepresentative, &nonPathRepresentative, target,
                                SPQR_INVALID_ARC));

   assert(
      target != pathRepresentative && target != nonPathRepresentative && pathRepresentative != nonPathRepresentative);
   assert(SPQRarcIsValid(pathRepresentative) && SPQRarcIsValid(nonPathRepresentative) && SPQRarcIsValid(target));
   assert(getNumMemberArcs(dec, member) == 3);

   //Create nodes
   spqr_node a = SPQR_INVALID_NODE;
   spqr_node b = SPQR_INVALID_NODE;
   spqr_node c = SPQR_INVALID_NODE;
   SCIP_CALL(createNode(dec, &a));
   SCIP_CALL(createNode(dec, &b));
   SCIP_CALL(createNode(dec, &c));

   // a -- b
   //  \  /
   //   c

   //Set arc nodes
   //Set target from b to c,
   SCIP_Bool targetReversed = arcIsReversedNonRigid(dec, target);
   setArcHeadAndTail(dec, target, c, b);

   MemberPathType pathType = newCol->reducedMembers[reducedMember].pathType;
   assert(pathType == INTO_HEAD || pathType == OUT_HEAD);
   if( arcIsReversedNonRigid(dec, pathRepresentative) == targetReversed )
   {
      setArcHeadAndTail(dec, pathRepresentative, a, c);
   } else
   {
      setArcHeadAndTail(dec, pathRepresentative, c, a);
   }
   if( arcIsReversedNonRigid(dec, nonPathRepresentative) == targetReversed )
   {
      setArcHeadAndTail(dec, nonPathRepresentative, b, a);
   } else
   {
      setArcHeadAndTail(dec, nonPathRepresentative, a, b);
   }
   //setup signed union find; all arcs are placed are not reversed. We pick an arbitrary arc as 'root' arc for this skeleton
   arcSetReversed(dec, target, FALSE);
   arcSetReversed(dec, pathRepresentative, FALSE);
   arcSetReversed(dec, nonPathRepresentative, FALSE);
   arcSetRepresentative(dec, target, SPQR_INVALID_ARC);
   arcSetRepresentative(dec, pathRepresentative, target);
   arcSetRepresentative(dec, nonPathRepresentative, target);
   *representativeArc = target;

   if( pathType == INTO_HEAD )
   {
      setTerminalTail(newColInfo, a);
   } else
   {
      setTerminalHead(newColInfo, a);
   }

   *mergedMember = member;

   return SCIP_OKAY;
}

static SCIP_RETCODE transformAndMergeParallel(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   reduced_member_id current,
   reduced_member_id next,
   spqr_member nextMember,
   SCIP_Bool nextIsParent,
   spqr_arc* representativeArc,
   spqr_member* mergedMember
)
{
   //Split off edges not in current subtree to form one parallel (or one edge)
   spqr_member childParallel = INVALID_REDUCED_MEMBER;
   spqr_arc source = newCol->reducedMembers[next].pathSourceArc;
   spqr_arc target = newCol->reducedMembers[next].pathTargetArc;
   if( getNumMemberArcs(dec, nextMember) > 3 )
   {
      SCIP_CALL(splitParallel(dec, nextMember, source, target, &childParallel));
      nextMember = childParallel;
      newCol->reducedMembers[next].member = childParallel;
   }
   assert(getNumMemberArcs(dec, nextMember) == 3);

   spqr_node sourceHead = SPQR_INVALID_NODE;
   spqr_node sourceTail = SPQR_INVALID_NODE;
   SCIP_CALL(createNode(dec, &sourceHead));
   SCIP_CALL(createNode(dec, &sourceTail));

   //set edge nodes and arc union-find data
   {
      spqr_arc firstArc = getFirstMemberArc(dec, nextMember);
      spqr_arc arc = firstArc;

      SCIP_Bool sourceReversed = arcIsReversedNonRigid(dec, source);
      do
      {
         if( arcIsReversedNonRigid(dec, arc) == sourceReversed )
         {
            setArcHeadAndTail(dec, arc, sourceHead, sourceTail);
         } else
         {
            setArcHeadAndTail(dec, arc, sourceTail, sourceHead);
         }
         arcSetRepresentative(dec, arc, source);
         arcSetReversed(dec, arc, FALSE);

         arc = getNextMemberArc(dec, arc);
      } while( arc != firstArc );

      arcSetRepresentative(dec, source, SPQR_INVALID_ARC);
   }

   //fix arc orientations of members; we cannot reflect for parallels.
   *representativeArc = mergeArcSigns(dec, *representativeArc, source, FALSE);

   spqr_member newMergedMember = SPQR_INVALID_MEMBER;
   if( nextIsParent )
   {
      SCIP_CALL( mergeGivenMemberIntoParent(dec, *mergedMember, nextMember,
                                 source, newCol->reducedMembers[current].pathTargetArc, TRUE, &newMergedMember) );
   } else
   {
      SCIP_CALL( mergeGivenMemberIntoParent(dec, nextMember, *mergedMember,
                                 newCol->reducedMembers[current].pathTargetArc, source, TRUE, &newMergedMember) );
   }
   *mergedMember = newMergedMember;

   return SCIP_OKAY;
}

static SCIP_RETCODE transformAndMergeSeries(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   reduced_member_id current,
   reduced_member_id next,
   spqr_member nextMember,
   SCIP_Bool nextIsParent,
   spqr_arc* representativeArc,
   spqr_member* mergedMember,
   NewColInformation* info
)
{

   SPQRColReducedMember* redMem = &newCol->reducedMembers[next];
   spqr_arc source = redMem->pathSourceArc;
   spqr_arc target = redMem->pathTargetArc;

   spqr_arc pathRepresentative = SPQR_INVALID_ARC;
   spqr_arc nonPathRepresentative = SPQR_INVALID_ARC;
   SCIP_CALL(
      splitSeriesMerging(dec, newCol, redMem, nextMember, &pathRepresentative, &nonPathRepresentative, source, target));
   //After splitting there is the following possibilities for nodes a-d:
   //(a)-source-(b)-path-(c)-target-(d)-nonpath-(a)
   //(a)-source-(b)-path-(c)-target-(d==a)
   //(a)-source-(b)=(c)-target-(d)-nonpath-(a)
   //(a)-source-(b)-path-(c)=(d) -nonpath-(a)
   //Note that the given arc is always between the same nodes
   assert(getNumMemberArcs(dec, nextMember) == 3 || getNumMemberArcs(dec, nextMember) == 4);
   assert(pathRepresentative != source && nonPathRepresentative != source &&
          ( SPQRarcIsInvalid(target) || ( target != pathRepresentative && target != nonPathRepresentative )));
   spqr_node a = SPQR_INVALID_NODE;
   spqr_node b = SPQR_INVALID_NODE;
   spqr_node c = SPQR_INVALID_NODE;
   spqr_node d = SPQR_INVALID_NODE;
   SCIP_CALL(createNode(dec, &a));
   SCIP_CALL(createNode(dec, &b));
   if( SPQRarcIsValid(pathRepresentative))
   {
      SCIP_CALL(createNode(dec, &c));
   } else
   {
      c = b;
   }
   SCIP_Bool hasNonPath = SPQRarcIsValid(nonPathRepresentative);
   SCIP_Bool hasTarget = SPQRarcIsValid(target);
   if( hasNonPath && hasTarget )
   {
      SCIP_CALL(createNode(dec, &d));
   } else
   {
      if( hasNonPath )
      {
         d = c;
      } else
      {
         d = a;
      }
   }

   SCIP_Bool sourceReversed = arcIsReversedNonRigid(dec, source);
   SCIP_Bool pathStartInHead = isHead(newCol->reducedMembers[current].pathType);
   if( pathStartInHead )
   {
      setArcHeadAndTail(dec, source, b, a);
   } else
   {
      setArcHeadAndTail(dec, source, a, b);
   }
   if( SPQRarcIsValid(pathRepresentative))
   {
      if(( arcIsReversedNonRigid(dec, pathRepresentative) == sourceReversed ) == pathStartInHead )
      {
         setArcHeadAndTail(dec, pathRepresentative, c, b);
      } else
      {
         setArcHeadAndTail(dec, pathRepresentative, b, c);
      }
      arcSetReversed(dec, pathRepresentative, FALSE);
      arcSetRepresentative(dec, pathRepresentative, source);
   }
   if( hasTarget )
   {
      if(( arcIsReversedNonRigid(dec, target) == sourceReversed ) == pathStartInHead )
      {
         setArcHeadAndTail(dec, target, d, c);
      } else
      {
         setArcHeadAndTail(dec, target, c, d);
      }
      arcSetReversed(dec, target, FALSE);
      arcSetRepresentative(dec, target, source);
   }
   if( hasNonPath )
   {
      if(( arcIsReversedNonRigid(dec, nonPathRepresentative) == sourceReversed ) == pathStartInHead )
      {
         setArcHeadAndTail(dec, nonPathRepresentative, a, d);
      } else
      {
         setArcHeadAndTail(dec, nonPathRepresentative, d, a);
      }
      arcSetReversed(dec, nonPathRepresentative, FALSE);
      arcSetRepresentative(dec, nonPathRepresentative, source);
   }

   arcSetReversed(dec, source, FALSE);
   arcSetRepresentative(dec, source, SPQR_INVALID_ARC);

   //fix arc orientations of members; we cannot reflect for series

   spqr_member newMergedMember = SPQR_INVALID_MEMBER;
   if( nextIsParent )
   {
      SCIP_CALL( mergeGivenMemberIntoParent(dec, *mergedMember, nextMember,
                                 source, newCol->reducedMembers[current].pathTargetArc, TRUE, &newMergedMember) );
   } else
   {
      SCIP_CALL( mergeGivenMemberIntoParent(dec, nextMember, *mergedMember,
                                 newCol->reducedMembers[current].pathTargetArc, source, TRUE, &newMergedMember) );
   }
   *mergedMember = newMergedMember;

   *representativeArc = mergeArcSigns(dec, *representativeArc, source, FALSE);
   if( !hasTarget )
   {
      //We are in the last node; finish the path
      setTerminalReversed(info, FALSE);
      if( isInto(newCol->reducedMembers[current].pathType))
      {
         setTerminalHead(info, c);
      } else
      {
         setTerminalTail(info, c);
      }
      setTerminalMember(info, *mergedMember);
      setTerminalRepresentative(info, *representativeArc);
   }
   return SCIP_OKAY;
}

static SCIP_RETCODE transformAndMergeRigid(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   reduced_member_id current,
   reduced_member_id next,
   spqr_member nextMember,
   SCIP_Bool nextIsParent,
   spqr_arc* representativeArc,
   spqr_member* mergedMember,
   NewColInformation* info
)
{
   SPQRColReducedMember* redMem = &newCol->reducedMembers[next];
   spqr_arc source = redMem->pathSourceArc;
   spqr_arc sourceRepresentative = findArcSign(dec, source).representative;


   spqr_member newMergedMember = SPQR_INVALID_MEMBER;

   if( nextIsParent )
   {
      SCIP_CALL( mergeGivenMemberIntoParent(dec, *mergedMember, nextMember,
                                 source, newCol->reducedMembers[current].pathTargetArc, !redMem->reverseArcs,
                                 &newMergedMember) );
   } else
   {
      SCIP_CALL( mergeGivenMemberIntoParent(dec, nextMember, *mergedMember,
                                 newCol->reducedMembers[current].pathTargetArc, source, !redMem->reverseArcs,
                                 &newMergedMember) );
   }

   *mergedMember = newMergedMember;

   *representativeArc = mergeArcSigns(dec, *representativeArc, sourceRepresentative, redMem->reverseArcs);

   if( SPQRarcIsInvalid(redMem->pathTargetArc))
   {
      //We are in the last node; finish the path
      setTerminalReversed(info, FALSE);
      if( isInto(newCol->reducedMembers[current].pathType))
      {
         if( redMem->reverseArcs )
         {
            setTerminalHead(info, redMem->rigidPathStart);
         } else
         {
            setTerminalHead(info, redMem->rigidPathEnd);
         }
      } else
      {
         if( redMem->reverseArcs )
         {
            setTerminalTail(info, redMem->rigidPathEnd);
         } else
         {
            setTerminalTail(info, redMem->rigidPathStart);
         }
      }
      setTerminalMember(info, *mergedMember);
      setTerminalRepresentative(info, *representativeArc);
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE transformAndMerge(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   reduced_member_id current,
   reduced_member_id next,
   spqr_arc* representativeArc,
   spqr_member* mergedMember,
   SCIP_Bool nextIsParent,
   NewColInformation* info
)
{
   spqr_member nextMember = newCol->reducedMembers[next].member;
   switch( getMemberType(dec, nextMember))
   {
      case SPQR_MEMBERTYPE_RIGID:
      {
         SCIP_CALL(transformAndMergeRigid(dec, newCol, current, next, nextMember, nextIsParent,
                                          representativeArc, mergedMember, info));
         break;
      }
      case SPQR_MEMBERTYPE_PARALLEL:
      {
         SCIP_CALL(transformAndMergeParallel(dec, newCol, current, next, nextMember, nextIsParent,
                                             representativeArc, mergedMember));
         break;
      }
      case SPQR_MEMBERTYPE_SERIES:
      {
         SCIP_CALL(transformAndMergeSeries(dec, newCol, current, next, nextMember, nextIsParent,
                                           representativeArc, mergedMember, info));
         break;
      }
      case SPQR_MEMBERTYPE_LOOP:
      case SPQR_MEMBERTYPE_UNASSIGNED:
      {
         assert(FALSE);
         break;
      }
   }
   return SCIP_OKAY;
}

static SCIP_RETCODE transformPath(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   SPQRColReducedComponent* component,
   NewColInformation* newColInfo
)
{
   //Realize first member
   reduced_member_id firstPathMember = component->pathEndMembers[0];

   spqr_arc representativeArc = SPQR_INVALID_ARC;
   spqr_member mergedMember = SPQR_INVALID_MEMBER;
   SCIP_CALL(transformFirstPathMember(dec, newCol, firstPathMember, newColInfo, &representativeArc, &mergedMember));
   //Iteratively call function which realizes next member and merges them together.
   reduced_member_id current = firstPathMember;
   reduced_member_id next = newCol->reducedMembers[current].nextPathMember;
   SCIP_Bool nextIsParent = newCol->reducedMembers[current].nextPathMemberIsParent;
   while( reducedMemberIsValid(next))
   {
      SCIP_CALL(
         transformAndMerge(dec, newCol, current, next, &representativeArc, &mergedMember, nextIsParent, newColInfo));
      current = next;
      next = newCol->reducedMembers[current].nextPathMember;
      nextIsParent = newCol->reducedMembers[current].nextPathMemberIsParent;
   }
   return SCIP_OKAY;
}

static SCIP_RETCODE columnTransformSingleParallel(
   SCIP_NETCOLADD* newCol,
   reduced_member_id reducedMemberId,
   spqr_member member,
   NewColInformation* newColInfo
)
{
   SPQRColReducedMember* reducedMember = &newCol->reducedMembers[reducedMemberId];
   assert(pathArcIsValid(reducedMember->firstPathArc) && reducedMember->numPathArcs == 1);
   //The new arc can be placed in parallel; just add it to this member
   setTerminalReversed(newColInfo, reducedMember->pathBackwards);
   setTerminalMember(newColInfo, member);
   return SCIP_OKAY;
}

static SCIP_RETCODE columnTransformSingleSeries(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   reduced_member_id reducedMemberId,
   spqr_member member,
   NewColInformation* newColInfo
)
{

   if( getNumMemberArcs(dec, member) == 1 )
   {
      newColInfo->member = member;
      newColInfo->reversed = newCol->arcInPathReversed[getFirstMemberArc(dec, member)];
      return SCIP_OKAY;
   }
   //Isolated single cycle
   spqr_member loopMember;
   SPQRColReducedMember* reducedMember = &newCol->reducedMembers[reducedMemberId];
   SCIP_CALL(splitSeries(dec, newCol, reducedMember, member, &loopMember, newColInfo));

   return SCIP_OKAY;
}

static SCIP_RETCODE columnTransformSingleRigid(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   reduced_member_id reducedMemberId,
   spqr_member member,
   NewColInformation* newColInfo
)
{
   assert(dec);
   assert(newCol);
   assert(newColInfo);
   assert(reducedMemberIsValid(reducedMemberId));
   //The path is already computed, so we can simply take the start and end nodes.
   //However, there is one exception, which is that an arc connecting these two nodes already exists in the member
   //If so, we create a new parallel member with the new arc and this member, unless the existing arc already points
   //to a parallel member
   SPQRColReducedMember* reducedMember = &newCol->reducedMembers[reducedMemberId];
   assert(SPQRnodeIsValid(reducedMember->rigidPathStart) && SPQRnodeIsValid(reducedMember->rigidPathEnd));
   {
      spqr_arc existingArcWithPath = SPQR_INVALID_ARC;
      spqr_arc firstArc = getFirstNodeArc(dec, reducedMember->rigidPathStart);
      spqr_arc arc = firstArc;
      SCIP_Bool pathInSameDirection = FALSE;
      do
      {
         spqr_node head = findArcHead(dec, arc);
         spqr_node tail = findArcTail(dec, arc);
         spqr_node other = head == reducedMember->rigidPathStart ? tail : head;
         if( other == reducedMember->rigidPathEnd )
         {
            existingArcWithPath = arc;
            pathInSameDirection = ( head == other ) != findArcSign(dec, existingArcWithPath).reversed;
            break;
         }
         arc = getNextNodeArc(dec, arc, reducedMember->rigidPathStart);
      } while( arc != firstArc );
      if( SPQRarcIsValid(existingArcWithPath))
      {
         SCIP_Bool isParent = FALSE;
         spqr_member adjacentMember = arcIsMarker(dec, existingArcWithPath) ?
            findArcChildMember(dec,existingArcWithPath) : SPQR_INVALID_MEMBER;
         if( existingArcWithPath == markerToParent(dec, member))
         {
            adjacentMember = findMemberParent(dec, member);
            isParent = TRUE;
         }
         if( SPQRmemberIsValid(adjacentMember) && getMemberType(dec, adjacentMember) == SPQR_MEMBERTYPE_PARALLEL )
         {
            spqr_arc parallelMarker = isParent ? markerOfParent(dec, member) : markerToParent(dec, adjacentMember);
            SCIP_Bool markerReversed = arcIsReversedNonRigid(dec, parallelMarker);
            setTerminalMember(newColInfo, adjacentMember);
            setTerminalReversed(newColInfo, markerReversed == pathInSameDirection);
         } else
         {
            //create a new parallel and move the edge there
            //This is a bit painful, because we cannot actually remove edges because of the union-find data structure
            //So what we do instead, is convert the current edge to a marker edge, and 'duplicate'
            //it in the new parallel member, and add the new marker there too, manually
            spqr_member adjacentParallel = SPQR_INVALID_MEMBER;
            SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &adjacentParallel));
            //'duplicate' a new arc in the parallel to be the current arc
            spqr_arc duplicate = SPQR_INVALID_ARC;
            spqr_element element = arcGetElement(dec, existingArcWithPath);
            if( element != MARKER_COLUMN_ELEMENT && element != MARKER_ROW_ELEMENT)
            {
               if( SPQRelementIsColumn(element))
               {
                  SCIP_CALL(createColumnArc(dec, adjacentParallel, &duplicate, SPQRelementToColumn(element), FALSE));
               } else
               {
                  SCIP_CALL(createRowArc(dec, adjacentParallel, &duplicate, SPQRelementToRow(element), FALSE));
               }
            } else if( isParent )
            {
               SCIP_CALL(createParentMarker(dec, adjacentParallel, arcIsTree(dec, existingArcWithPath), adjacentMember,
                                            markerOfParent(dec, member), &duplicate, FALSE));
            } else
            {
               SCIP_CALL(createChildMarker(dec, adjacentParallel, adjacentMember, arcIsTree(dec, existingArcWithPath),
                                           &duplicate, FALSE));
               dec->members[adjacentMember].parentMember = adjacentParallel;
               dec->members[adjacentMember].markerOfParent = duplicate;
            }
            //Create the other marker edge
            spqr_arc parallelMarker = SPQR_INVALID_ARC;
            if( isParent )
            {
               SCIP_CALL(createChildMarker(dec, adjacentParallel, member, !arcIsTree(dec, existingArcWithPath),
                                           &parallelMarker, FALSE));
            } else
            {
               SCIP_CALL(createParentMarker(dec, adjacentParallel, !arcIsTree(dec, existingArcWithPath),
                                            member, existingArcWithPath, &parallelMarker, FALSE));
            }

            //Change the existing edge to a marker
            if( isParent )
            {
               assert(markerToParent(dec, member) == existingArcWithPath);
               dec->arcs[markerOfParent(dec, member)].childMember = adjacentParallel;
               dec->members[member].parentMember = adjacentParallel;
               dec->members[member].markerToParent = existingArcWithPath;
               dec->members[member].markerOfParent = parallelMarker;
               dec->arcs[existingArcWithPath].element = arcIsTree(dec, existingArcWithPath) ? MARKER_ROW_ELEMENT
                                                                                            : MARKER_COLUMN_ELEMENT;
               dec->arcs[existingArcWithPath].childMember = adjacentParallel;

            } else
            {
               dec->arcs[existingArcWithPath].element = arcIsTree(dec, existingArcWithPath) ? MARKER_ROW_ELEMENT
                                                                                            : MARKER_COLUMN_ELEMENT;
               dec->arcs[existingArcWithPath].childMember = adjacentParallel;
            }

            setTerminalMember(newColInfo, adjacentParallel);
            setTerminalReversed(newColInfo, !pathInSameDirection);
         }
         return SCIP_OKAY;
      }
   }

   setTerminalMember(newColInfo, member);
   setTerminalReversed(newColInfo, FALSE);
   setTerminalTail(newColInfo, reducedMember->rigidPathStart);
   setTerminalHead(newColInfo, reducedMember->rigidPathEnd);
   setTerminalRepresentative(newColInfo,
                             findArcSign(dec, newCol->pathArcs[reducedMember->firstPathArc].arc).representative);
   return SCIP_OKAY;
}

static SCIP_RETCODE transformComponent(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol,
   SPQRColReducedComponent* component,
   NewColInformation* newColInfo
)
{
   assert(dec);
   assert(newCol);
   assert(component);
   assert(newColInfo);

   if( newCol->reducedMembers[component->root].numChildren ==
       newCol->reducedMembers[component->root].numPropagatedChildren )
   {
      //No merging necessary, only a single component
      reduced_member_id reducedMember = component->root;
      assert(reducedMemberIsValid(reducedMember));
      spqr_member member = newCol->reducedMembers[reducedMember].member;
      SPQRMemberType type = getMemberType(dec, member);

      switch( type )
      {
         case SPQR_MEMBERTYPE_RIGID:
         {
            SCIP_CALL(columnTransformSingleRigid(dec, newCol, reducedMember, member, newColInfo));
            break;
         }
         case SPQR_MEMBERTYPE_PARALLEL:
         {
            SCIP_CALL(columnTransformSingleParallel(newCol, reducedMember, member, newColInfo));
            break;
         }
         case SPQR_MEMBERTYPE_LOOP:
         case SPQR_MEMBERTYPE_SERIES:
         {
            SCIP_CALL(columnTransformSingleSeries(dec, newCol, reducedMember, member, newColInfo));
            break;
         }
         case SPQR_MEMBERTYPE_UNASSIGNED:
         default:
         {
            assert(FALSE);
            break;
         }
      }
      return SCIP_OKAY;
   }
   // Otherwise, the reduced members form a path which can be merged into a single component of type R
   SCIP_CALL(transformPath(dec, newCol, component, newColInfo));

   return SCIP_OKAY;
}
static SCIP_Bool SCIPnetcoladdRemainsNetwork(const SCIP_NETCOLADD* newCol)
{
   return newCol->remainsNetwork;
}

static SCIP_RETCODE SCIPnetcoladdAdd(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETCOLADD* newCol
)
{
   assert(dec);
   assert(newCol);
   assert(SCIPnetcoladdRemainsNetwork(newCol));

   if( newCol->numReducedComponents == 0 )
   {
      spqr_member member;
      SCIP_CALL(createStandaloneSeries(dec, newCol->newRowArcs, newCol->newRowArcReversed,
                                       newCol->numNewRowArcs, newCol->newColIndex, &member));
   } else if( newCol->numReducedComponents == 1 )
   {
      NewColInformation information = emptyNewColInformation();
      SCIP_CALL(transformComponent(dec, newCol, &newCol->reducedComponents[0], &information));
      assert(memberIsRepresentative(dec, information.member));
      if( newCol->numNewRowArcs == 0 )
      {
         spqr_arc colArc = SPQR_INVALID_ARC;
         SCIP_CALL(createColumnArc(dec, information.member, &colArc, newCol->newColIndex, information.reversed));
         if( SPQRnodeIsValid(information.head))
         {
            assert(SPQRnodeIsValid(information.tail));
            assert(SPQRarcIsValid(information.representative));
            setArcHeadAndTail(dec, colArc, findNode(dec, information.head), findNode(dec, information.tail));
            arcSetRepresentative(dec, colArc, information.representative);
            arcSetReversed(dec, colArc, information.reversed != arcIsReversedNonRigid(dec, information.representative));
         }
      } else
      {
         spqr_member newSeries = SPQR_INVALID_MEMBER;
         SCIP_CALL(createConnectedSeries(dec, newCol->newRowArcs, newCol->newRowArcReversed, newCol->numNewRowArcs,
                                         newCol->newColIndex, &newSeries));
         spqr_arc markerArc = SPQR_INVALID_ARC;
         spqr_arc ignore = SPQR_INVALID_ARC;
         SCIP_CALL(createMarkerPairWithReferences(dec, information.member, newSeries, FALSE, information.reversed, TRUE,
                                                  &markerArc, &ignore));
         if( SPQRnodeIsValid(information.head))
         {
            assert(SPQRnodeIsValid(information.tail));
            assert(SPQRarcIsValid(information.representative));
            setArcHeadAndTail(dec, markerArc, findNode(dec, information.head), findNode(dec, information.tail));
            arcSetRepresentative(dec, markerArc, information.representative);
            arcSetReversed(dec, markerArc,
                           information.reversed != arcIsReversedNonRigid(dec, information.representative));
         }
      }
      if( getMemberType(dec, information.member) == SPQR_MEMBERTYPE_LOOP )
      {
         assert(getNumMemberArcs(dec, information.member) == 2 || getNumMemberArcs(dec, information.member) == 3);
         if( getNumMemberArcs(dec, information.member) == 3 )
         {
            changeLoopToParallel(dec, information.member);
         }
      }
   } else
   {
#ifndef NDEBUG
      int numDecComponentsBefore = numConnectedComponents(dec);
#endif
      spqr_member newSeries = SPQR_INVALID_MEMBER;
      SCIP_CALL(createConnectedSeries(dec, newCol->newRowArcs, newCol->newRowArcReversed,
                                      newCol->numNewRowArcs, newCol->newColIndex, &newSeries));
      for( int i = 0; i < newCol->numReducedComponents; ++i )
      {
         NewColInformation information = emptyNewColInformation();
         SCIP_CALL(transformComponent(dec, newCol, &newCol->reducedComponents[i], &information));
         if( getMemberType(dec, information.member) == SPQR_MEMBERTYPE_LOOP )
         {
            assert(getNumMemberArcs(dec, information.member) == 1);
            spqr_arc arc = getFirstMemberArc(dec, information.member);
            assert(newCol->arcInPath[arc]);
            moveArcToNewMember(dec, arc, information.member, newSeries);
            arcSetReversed(dec, arc, !newCol->arcInPathReversed[arc]);
            dec->members[information.member].type = SPQR_MEMBERTYPE_UNASSIGNED;
         } else
         {
            reorderComponent(dec,
                             information.member);//reorder the subtree so that the newly series member is a parent
            spqr_arc markerArc = SPQR_INVALID_ARC;
            spqr_arc ignore = SPQR_INVALID_ARC;
            SCIP_CALL(createMarkerPairWithReferences(dec, newSeries, information.member, TRUE, information.reversed,
                                                     TRUE, &ignore, &markerArc));
            if( SPQRnodeIsValid(information.head))
            {
               assert(SPQRnodeIsValid(information.tail));
               assert(SPQRarcIsValid(information.representative));
               setArcHeadAndTail(dec, markerArc, findNode(dec, information.head), findNode(dec, information.tail));
               arcSetRepresentative(dec, markerArc, information.representative);
               arcSetReversed(dec, markerArc,
                              information.reversed == arcIsReversedNonRigid(dec, information.representative));
            }
         }
      }
      decreaseNumConnectedComponents(dec, newCol->numReducedComponents - 1);
      assert(numConnectedComponents(dec) == ( numDecComponentsBefore - newCol->numReducedComponents + 1 ));
   }
   return SCIP_OKAY;
}




static int minValue(
   int a,
   int b
)
{
   return a < b ? a : b;
}

typedef int cut_arc_id;
#define INVALID_CUT_ARC (-1)

static SCIP_Bool cutArcIsInvalid(const cut_arc_id arc)
{
   return arc < 0;
}

static SCIP_Bool cutArcIsValid(const cut_arc_id arc)
{
   return !cutArcIsInvalid(arc);
}

typedef struct
{//TODO:test if memory overhead of pointers is worth it?
   spqr_arc arc;
   spqr_node arcHead;
   spqr_node arcTail;
   cut_arc_id nextMember;
   cut_arc_id nextOverall;
   SCIP_Bool arcReversed;
} CutArcListNode;


typedef enum
{
   TYPE_UNDETERMINED = 0,
   TYPE_PROPAGATED = 1,
   TYPE_MERGED = 2,
   TYPE_NOT_NETWORK = 3
} RowReducedMemberType;


typedef struct
{
   int low;
   int discoveryTime;
} ArticulationNodeInformation;

//We allocate the callstacks of recursive algorithms (usually DFS, bounded by some linear number of calls)
//If one does not do this, we overflow the stack for large matrices/graphs through the number of recursive function calls
//Then, we can write the recursive algorithms as while loops and allocate the function call data on the heap, preventing
//Stack overflows
typedef struct
{
   spqr_node node;
   spqr_arc nodeArc;
} DFSCallData;

typedef struct
{
   children_idx currentChild;
   reduced_member_id id;
} MergeTreeCallData;

typedef struct
{
   spqr_node node;
   spqr_arc arc;
} ColorDFSCallData;

typedef struct
{
   spqr_arc arc;
   spqr_node node;
   spqr_node parent;
   SCIP_Bool isAP;
} ArticulationPointCallStack;

typedef enum
{
   UNCOLORED = 0,
   COLOR_SOURCE = 1,
   COLOR_SINK = 2
} COLOR_STATUS;

typedef struct
{
   spqr_member member;
   spqr_member rootMember;
   int depth;
   RowReducedMemberType type;
   reduced_member_id parent;

   children_idx firstChild;
   children_idx numChildren;
   children_idx numPropagatedChildren;

   cut_arc_id firstCutArc;
   int numCutArcs;

   //For non-rigid members
   spqr_arc splitArc;
   SCIP_Bool splitHead;    //Otherwise the tail of this arc is split
   SCIP_Bool otherIsSource;//Otherwise the other end node is part of the sink partition.
   // For non-rigid members this refers to splitArc, for rigid members it refers to articulation arc

   //For rigid members
   spqr_node otherNode;
   spqr_node splitNode;
   SCIP_Bool allHaveCommonNode;
   SCIP_Bool otherNodeSplit;
   SCIP_Bool willBeReversed;
   spqr_arc articulationArc;
   spqr_node coloredNode;//points to a colored node so that we can efficiently zero out the colors again.

} SPQRRowReducedMember;

typedef struct
{
   int rootDepth;
   reduced_member_id root;
} SPQRRowReducedComponent;

typedef struct
{
   SCIP_Bool remainsNetwork;

   SPQRRowReducedMember* reducedMembers;
   int memReducedMembers;
   int numReducedMembers;

   SPQRRowReducedComponent* reducedComponents;
   int memReducedComponents;
   int numReducedComponents;

   MemberInfo* memberInformation;
   int memMemberInformation;
   int numMemberInformation;

   reduced_member_id* childrenStorage;
   int memChildrenStorage;
   int numChildrenStorage;

   CutArcListNode* cutArcs;
   int memCutArcs;
   int numCutArcs;
   cut_arc_id firstOverallCutArc;

   spqr_row newRowIndex;

   spqr_col* newColumnArcs;
   SCIP_Bool* newColumnReversed;
   int memColumnArcs;
   int numColumnArcs;

   reduced_member_id* leafMembers;
   int numLeafMembers;
   int memLeafMembers;

   spqr_arc* decompositionColumnArcs;
   SCIP_Bool* decompositionColumnArcReversed;
   int memDecompositionColumnArcs;
   int numDecompositionColumnArcs;

   SCIP_Bool* isArcCut;
   SCIP_Bool* isArcCutReversed;
   int numIsArcCut;
   int memIsArcCut;

   COLOR_STATUS* nodeColors;
   int memNodeColors;

   spqr_node* articulationNodes;
   int numArticulationNodes;
   int memArticulationNodes;

   ArticulationNodeInformation* articulationNodeSearchInfo;
   int memNodeSearchInfo;

   int* crossingPathCount;
   int memCrossingPathCount;

   DFSCallData* intersectionDFSData;
   int memIntersectionDFSData;

   ColorDFSCallData* colorDFSData;
   int memColorDFSData;

   ArticulationPointCallStack* artDFSData;
   int memArtDFSData;

   CreateReducedMembersCallstack* createReducedMembersCallstack;
   int memCreateReducedMembersCallstack;

   int* intersectionPathDepth;
   int memIntersectionPathDepth;

   spqr_node* intersectionPathParent;
   int memIntersectionPathParent;

   MergeTreeCallData* mergeTreeCallData;
   int memMergeTreeCallData;

   COLOR_STATUS * temporaryColorArray;
   int memTemporaryColorArray;
} SCIP_NETROWADD;

typedef struct
{
   spqr_member member;
   spqr_node head;
   spqr_node tail;
   spqr_arc representative;
   SCIP_Bool reversed;
} NewRowInformation;

static NewRowInformation emptyNewRowInformation(void)
{
   NewRowInformation information;
   information.member = SPQR_INVALID_MEMBER;
   information.head = SPQR_INVALID_NODE;
   information.tail = SPQR_INVALID_NODE;
   information.representative = SPQR_INVALID_ARC;
   information.reversed = FALSE;
   return information;
}

/**
 * Saves the information of the current row and partitions it based on whether or not the given columns are
 * already part of the decomposition.
 */
static SCIP_RETCODE newRowUpdateRowInformation(
   const SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const spqr_row row,
   const spqr_col* columns,
   const double* columnValues,
   const int numColumns
)
{
   newRow->newRowIndex = row;

   newRow->numDecompositionColumnArcs = 0;
   newRow->numColumnArcs = 0;

   for( int i = 0; i < numColumns; ++i )
   {
      spqr_arc columnArc = getDecompositionColumnArc(dec, columns[i]);
      SCIP_Bool reversed = columnValues[i] < 0.0;
      if( SPQRarcIsValid(columnArc))
      {//If the arc is the current decomposition: save it in the array
         if( newRow->numDecompositionColumnArcs == newRow->memDecompositionColumnArcs )
         {
            int newNumArcs = newRow->memDecompositionColumnArcs == 0 ? 8 : 2 * newRow->memDecompositionColumnArcs;
            SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->decompositionColumnArcs,
                                                  newRow->memDecompositionColumnArcs, newNumArcs));
            SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->decompositionColumnArcReversed,
                                                  newRow->memDecompositionColumnArcs, newNumArcs));
            newRow->memDecompositionColumnArcs = newNumArcs;
         }
         newRow->decompositionColumnArcs[newRow->numDecompositionColumnArcs] = columnArc;
         newRow->decompositionColumnArcReversed[newRow->numDecompositionColumnArcs] = reversed;
         ++newRow->numDecompositionColumnArcs;
      } else
      {
         //Not in the decomposition: add it to the set of arcs which are newly added with this row.
         if( newRow->numColumnArcs == newRow->memColumnArcs )
         {
            int newNumArcs = newRow->memColumnArcs == 0 ? 8 : 2 * newRow->memColumnArcs;
            SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->newColumnArcs,
                                                  newRow->memColumnArcs, newNumArcs));
            SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->newColumnReversed,
                                                  newRow->memColumnArcs, newNumArcs));
            newRow->memColumnArcs = newNumArcs;
         }
         newRow->newColumnArcs[newRow->numColumnArcs] = columns[i];
         newRow->newColumnReversed[newRow->numColumnArcs] = reversed;
         newRow->numColumnArcs++;
      }
   }

   return SCIP_OKAY;
}

/**
 * Recursively creates reduced members from this member to the root of the decomposition tree.
 * @param dec
 * @param newRow
 * @param member
 * @return
 */
static reduced_member_id createRowReducedMembersToRoot(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const spqr_member firstMember
)
{
   assert(SPQRmemberIsValid(firstMember));

   CreateReducedMembersCallstack* callstack = newRow->createReducedMembersCallstack;
   callstack[0].member = firstMember;
   int callDepth = 0;

   while( callDepth >= 0 )
   {
      spqr_member member = callstack[callDepth].member;
      reduced_member_id reducedMember = newRow->memberInformation[member].reducedMember;

      SCIP_Bool reducedValid = reducedMemberIsValid(reducedMember);
      if( !reducedValid )
      {
         //reduced member was not yet created; we create it
         reducedMember = newRow->numReducedMembers;

         SPQRRowReducedMember* reducedMemberData = &newRow->reducedMembers[reducedMember];
         ++newRow->numReducedMembers;

         reducedMemberData->member = member;
         reducedMemberData->numChildren = 0;
         reducedMemberData->numCutArcs = 0;
         reducedMemberData->firstCutArc = INVALID_CUT_ARC;
         reducedMemberData->type = TYPE_UNDETERMINED;
         reducedMemberData->numPropagatedChildren = 0;
         reducedMemberData->articulationArc = SPQR_INVALID_ARC;
         reducedMemberData->splitNode = SPQR_INVALID_NODE;
         reducedMemberData->otherNode = SPQR_INVALID_NODE;
         reducedMemberData->splitArc = SPQR_INVALID_ARC;
         reducedMemberData->splitHead = FALSE;
         reducedMemberData->allHaveCommonNode = FALSE;
         reducedMemberData->otherNodeSplit = FALSE;
         reducedMemberData->willBeReversed = FALSE;
         reducedMemberData->coloredNode = SPQR_INVALID_NODE;

         newRow->memberInformation[member].reducedMember = reducedMember;
         assert(memberIsRepresentative(dec, member));
         spqr_member parentMember = findMemberParent(dec, member);

         if( SPQRmemberIsValid(parentMember))
         {
            //recursive call to parent member
            ++callDepth;
            assert(callDepth < newRow->memCreateReducedMembersCallstack);
            callstack[callDepth].member = parentMember;
            continue;

         } else
         {
            //we found a new reduced decomposition component

            reducedMemberData->parent = INVALID_REDUCED_MEMBER;
            reducedMemberData->depth = 0;
            reducedMemberData->rootMember = member;

            assert(newRow->numReducedComponents < newRow->memReducedComponents);
            newRow->reducedComponents[newRow->numReducedComponents].root = reducedMember;
            ++newRow->numReducedComponents;
         }
      }
      if( reducedValid )
      {
         assert(reducedMember < newRow->numReducedMembers);
         //Reduced member was already created in earlier call
         //update the depth of the root if appropriate
         reduced_member_id* depthMinimizer = &newRow->memberInformation[newRow->reducedMembers[reducedMember].rootMember].rootDepthMinimizer;
         if( reducedMemberIsInvalid(*depthMinimizer) ||
             newRow->reducedMembers[reducedMember].depth < newRow->reducedMembers[*depthMinimizer].depth )
         {
            *depthMinimizer = reducedMember;
         }
      }
      while( TRUE ) /*lint !e716*/
      {
         --callDepth;
         if( callDepth < 0 ) break;
         spqr_member parentMember = callstack[callDepth + 1].member;
         reduced_member_id parentReducedMember = newRow->memberInformation[parentMember].reducedMember;
         spqr_member currentMember = callstack[callDepth].member;
         reduced_member_id currentReducedMember = newRow->memberInformation[currentMember].reducedMember;

         SPQRRowReducedMember* parentReducedMemberData = &newRow->reducedMembers[parentReducedMember];
         SPQRRowReducedMember* reducedMemberData = &newRow->reducedMembers[currentReducedMember];

         reducedMemberData->parent = parentReducedMember;
         reducedMemberData->depth = parentReducedMemberData->depth + 1;
         reducedMemberData->rootMember = parentReducedMemberData->rootMember;

         newRow->reducedMembers[parentReducedMember].numChildren++;
      }
   }

   reduced_member_id returnedMember = newRow->memberInformation[callstack[0].member].reducedMember;
   return returnedMember;
}


/**
 * Construct a smaller sub tree of the decomposition on which the cut arcs lie.
 * @return
 */
static SCIP_RETCODE constructRowReducedDecomposition(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow
)
{
#ifndef NDEBUG
   for( int i = 0; i < newRow->memMemberInformation; ++i )
   {
      assert(reducedMemberIsInvalid(newRow->memberInformation[i].reducedMember));
   }
#endif

   newRow->numReducedComponents = 0;
   newRow->numReducedMembers = 0;
   if( newRow->numDecompositionColumnArcs == 0 )
   {//Early return in case the reduced decomposition will be empty
      return SCIP_OKAY;
   }
   assert(newRow->numReducedMembers == 0);
   assert(newRow->numReducedComponents == 0);

   int newSize = largestMemberID(dec);//Is this sufficient?
   if( newSize > newRow->memReducedMembers )
   {
      int updatedSize = maxValue(2 * newRow->memReducedMembers, newSize);
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->reducedMembers, newRow->memReducedMembers, updatedSize));
      newRow->memReducedMembers = updatedSize;
   }
   if( newSize > newRow->memMemberInformation )
   {
      int updatedSize = maxValue(2 * newRow->memMemberInformation, newSize);
      SCIP_CALL(
         SCIPreallocBlockMemoryArray(dec->env, &newRow->memberInformation, newRow->memMemberInformation, updatedSize));
      for( int i = newRow->memMemberInformation; i < updatedSize; ++i )
      {
         newRow->memberInformation[i].reducedMember = INVALID_REDUCED_MEMBER;
         newRow->memberInformation[i].rootDepthMinimizer = INVALID_REDUCED_MEMBER;
      }
      newRow->memMemberInformation = updatedSize;
   }

   int numComponents = numConnectedComponents(dec);
   if( numComponents > newRow->memReducedComponents )
   {
      int updatedSize = maxValue(2 * newRow->memReducedComponents, numComponents);
      SCIP_CALL(
         SCIPreallocBlockMemoryArray(dec->env, &newRow->reducedComponents, newRow->memReducedComponents, updatedSize));
      newRow->memReducedComponents = updatedSize;
   }

   int numMembers = getNumMembers(dec);
   if( newRow->memCreateReducedMembersCallstack < numMembers )
   {
      int updatedSize = maxValue(2 * newRow->memCreateReducedMembersCallstack, numMembers);
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->createReducedMembersCallstack,
                                            newRow->memCreateReducedMembersCallstack, updatedSize));
      newRow->memCreateReducedMembersCallstack = updatedSize;
   }

   //Create the reduced members (recursively)
   for( int i = 0; i < newRow->numDecompositionColumnArcs; ++i )
   {
      assert(i < newRow->memDecompositionColumnArcs);
      spqr_arc arc = newRow->decompositionColumnArcs[i];
      spqr_member arcMember = findArcMember(dec, arc);
      reduced_member_id reducedMember = createRowReducedMembersToRoot(dec, newRow, arcMember);
      reduced_member_id* depthMinimizer = &newRow->memberInformation[newRow->reducedMembers[reducedMember].rootMember].rootDepthMinimizer;
      if( reducedMemberIsInvalid(*depthMinimizer))
      {
         *depthMinimizer = reducedMember;
      }
   }

   //Set the reduced roots according to the root depth minimizers
   for( int i = 0; i < newRow->numReducedComponents; ++i )
   {
      SPQRRowReducedComponent* component = &newRow->reducedComponents[i];
      spqr_member rootMember = newRow->reducedMembers[component->root].member;
      reduced_member_id reducedMinimizer = newRow->memberInformation[rootMember].rootDepthMinimizer;
      component->rootDepth = newRow->reducedMembers[reducedMinimizer].depth;
      component->root = reducedMinimizer;

      //This simplifies code further down which does not need to be component-aware; just pretend that the reduced member is the new root.
      newRow->reducedMembers[component->root].parent = INVALID_REDUCED_MEMBER;
      assert(memberIsRepresentative(dec, rootMember));
   }

   //update the children array
   int numTotalChildren = 0;
   for( int i = 0; i < newRow->numReducedMembers; ++i )
   {
      SPQRRowReducedMember* reducedMember = &newRow->reducedMembers[i];
      reduced_member_id minimizer = newRow->memberInformation[reducedMember->rootMember].rootDepthMinimizer;
      if( reducedMember->depth >= newRow->reducedMembers[minimizer].depth )
      {
         reducedMember->firstChild = numTotalChildren;
         numTotalChildren += reducedMember->numChildren;
         reducedMember->numChildren = 0;
      }
   }

   if( newRow->memChildrenStorage < numTotalChildren )
   {
      int newMemSize = maxValue(newRow->memChildrenStorage * 2, numTotalChildren);
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->childrenStorage, newRow->memChildrenStorage,
                                            newMemSize));
      newRow->memChildrenStorage = newMemSize;
   }
   newRow->numChildrenStorage = numTotalChildren;

   //Fill up the children array`
   for( reduced_member_id reducedMember = 0; reducedMember < newRow->numReducedMembers; ++reducedMember )
   {
      SPQRRowReducedMember* reducedMemberData = &newRow->reducedMembers[reducedMember];
      if( reducedMemberData->depth <=
          newRow->reducedMembers[newRow->memberInformation[reducedMemberData->rootMember].rootDepthMinimizer].depth )
      {
         continue;
      }
      spqr_member parentMember = findMemberParent(dec, reducedMemberData->member);
      reduced_member_id parentReducedMember = SPQRmemberIsValid(parentMember)
                                              ? newRow->memberInformation[parentMember].reducedMember
                                              : INVALID_REDUCED_MEMBER;
      if( reducedMemberIsValid(parentReducedMember))
      {
         SPQRRowReducedMember* parentReducedMemberData = &newRow->reducedMembers[parentReducedMember];
         newRow->childrenStorage[parentReducedMemberData->firstChild +
                                 parentReducedMemberData->numChildren] = reducedMember;
         ++parentReducedMemberData->numChildren;
      }
   }

   //Clean up the root depth minimizers.
   for( int i = 0; i < newRow->numReducedMembers; ++i )
   {
      SPQRRowReducedMember* reducedMember = &newRow->reducedMembers[i];
      assert(reducedMember);
      spqr_member rootMember = reducedMember->rootMember;
      assert(rootMember >= 0);
      assert(rootMember < dec->memMembers);
      newRow->memberInformation[rootMember].rootDepthMinimizer = INVALID_REDUCED_MEMBER;
   }

   return SCIP_OKAY;
}


/**
 * Marks an arc as 'cut'. This implies that its cycle in the decomposition must be elongated
 * @param newRow
 * @param arc
 * @param reducedMember
 */
static void createCutArc(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const spqr_arc arc,
   const reduced_member_id reducedMember,
   SCIP_Bool reversed
)
{
   cut_arc_id cut_arc = newRow->numCutArcs;

   CutArcListNode* listNode = &newRow->cutArcs[cut_arc];
   listNode->arc = arc;

   listNode->nextMember = newRow->reducedMembers[reducedMember].firstCutArc;
   newRow->reducedMembers[reducedMember].firstCutArc = cut_arc;

   listNode->nextOverall = newRow->firstOverallCutArc;
   newRow->firstOverallCutArc = cut_arc;

   newRow->numCutArcs++;
   newRow->reducedMembers[reducedMember].numCutArcs++;
   assert(newRow->numCutArcs <= newRow->memCutArcs);

   assert(arc < newRow->memIsArcCut);
   newRow->isArcCut[arc] = TRUE;
   newRow->isArcCutReversed[arc] = reversed;

   assert(memberIsRepresentative(dec, newRow->reducedMembers[reducedMember].member));
   if( getMemberType(dec, newRow->reducedMembers[reducedMember].member) == SPQR_MEMBERTYPE_RIGID )
   {

      listNode->arcHead = findEffectiveArcHead(dec, arc);
      listNode->arcTail = findEffectiveArcTail(dec, arc);
      if( reversed )
      {
         SCIPswapInts(&listNode->arcHead, &listNode->arcTail);
      }
      assert(SPQRnodeIsValid(listNode->arcHead) && SPQRnodeIsValid(listNode->arcTail));
   } else
   {
      listNode->arcHead = SPQR_INVALID_NODE;
      listNode->arcTail = SPQR_INVALID_NODE;
   }

   listNode->arcReversed = reversed;
}

/**
 * Creates all cut arcs within the decomposition for the new row.
 * Note this preallocates memory for cut arcs which may be created by propagation.
 */
static SCIP_RETCODE createReducedDecompositionCutArcs(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow
)
{
   //Allocate memory for cut arcs
   spqr_arc maxArcID = largestArcID(dec);
   if( maxArcID > newRow->memIsArcCut )
   {
      int newSize = maxValue(maxArcID, 2 * newRow->memIsArcCut);
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->isArcCut, newRow->memIsArcCut, newSize));
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->isArcCutReversed, newRow->memIsArcCut, newSize));
      for( int i = newRow->memIsArcCut; i < newSize; ++i )
      {
         newRow->isArcCut[i] = FALSE;
         newRow->isArcCutReversed[i] = FALSE;
      }
      newRow->memIsArcCut = newSize;
   }
#ifndef NDEBUG
   for( int i = 0; i < newRow->memIsArcCut; ++i )
   {
      assert(!newRow->isArcCut[i]);
      assert(!newRow->isArcCutReversed[i]);
   }
#endif

   int numNeededArcs = newRow->numDecompositionColumnArcs * 4;
   if( numNeededArcs > newRow->memCutArcs )
   {
      int newSize = maxValue(newRow->memCutArcs * 2, numNeededArcs);
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->cutArcs, newRow->memCutArcs, newSize));
      newRow->memCutArcs = newSize;
   }
   newRow->numCutArcs = 0;
   newRow->firstOverallCutArc = INVALID_CUT_ARC;
   for( int i = 0; i < newRow->numDecompositionColumnArcs; ++i )
   {
      spqr_arc arc = newRow->decompositionColumnArcs[i];
      spqr_member member = findArcMember(dec, arc);
      reduced_member_id reduced_member = newRow->memberInformation[member].reducedMember;
      assert(reducedMemberIsValid(reduced_member));
      createCutArc(dec, newRow, arc, reduced_member, newRow->decompositionColumnArcReversed[i]);
   }

   return SCIP_OKAY;
}

/**
 * Determines the members of the reduced decomposition which are leafs.
 * This is used in propagation to ensure propagation is only checked for components which have at most one neighbour
 * which is not propagated.
 */
static SCIP_RETCODE determineLeafReducedMembers(
   const SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow
)
{
   if( newRow->numDecompositionColumnArcs > newRow->memLeafMembers )
   {
      int updatedSize = maxValue(newRow->numDecompositionColumnArcs, 2 * newRow->memLeafMembers);
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->leafMembers, newRow->memLeafMembers, updatedSize));
      newRow->memLeafMembers = updatedSize;
   }
   newRow->numLeafMembers = 0;

   for( reduced_member_id reducedMember = 0; reducedMember < newRow->numReducedMembers; ++reducedMember )
   {
      if( newRow->reducedMembers[reducedMember].numChildren == 0 )
      {
         newRow->leafMembers[newRow->numLeafMembers] = reducedMember;
         ++newRow->numLeafMembers;
      }
   }
   return SCIP_OKAY;
}

/**
 * Preallocates memory arrays necessary for searching rigid components.
 */
static SCIP_RETCODE allocateRigidSearchMemory(
   const SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow
)
{
   int totalNumNodes = getNumNodes(dec);
   int maxNumNodes = 2 * dec->numArcs;
   if( maxNumNodes > newRow->memNodeColors )
   {
      int newSize = maxValue(2 * newRow->memNodeColors, maxNumNodes);
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->nodeColors, newRow->memNodeColors, newSize));
      for( int i = newRow->memNodeColors; i < newSize; ++i )
      {
         newRow->nodeColors[i] = UNCOLORED;
      }
      newRow->memNodeColors = newSize;
   }

   if( totalNumNodes > newRow->memArticulationNodes )
   {
      int newSize = maxValue(2 * newRow->memArticulationNodes, totalNumNodes);
      SCIP_CALL(
         SCIPreallocBlockMemoryArray(dec->env, &newRow->articulationNodes, newRow->memArticulationNodes, newSize));
      newRow->memArticulationNodes = newSize;
   }
   if( totalNumNodes > newRow->memNodeSearchInfo )
   {
      int newSize = maxValue(2 * newRow->memNodeSearchInfo, totalNumNodes);
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->articulationNodeSearchInfo, newRow->memNodeSearchInfo,
                                            newSize));
      newRow->memNodeSearchInfo = newSize;
   }
   if( totalNumNodes > newRow->memCrossingPathCount )
   {
      int newSize = maxValue(2 * newRow->memCrossingPathCount, totalNumNodes);
      SCIP_CALL(
         SCIPreallocBlockMemoryArray(dec->env, &newRow->crossingPathCount, newRow->memCrossingPathCount, newSize));
      newRow->memCrossingPathCount = newSize;
   }
   if( totalNumNodes > newRow->memTemporaryColorArray ){
      int newSize = maxValue(2 * newRow->memTemporaryColorArray, totalNumNodes);
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->temporaryColorArray,
                                            newRow->memTemporaryColorArray, newSize));
      newRow->memTemporaryColorArray = newSize;
   }

   //TODO: see if tradeoff for performance bound by checking max # of nodes of rigid is worth it to reduce size
   //of the following allocations
   int largestID = largestNodeID(
      dec);//TODO: only update the stack sizes of the following when needed? The preallocation might be causing performance problems
   if( largestID > newRow->memIntersectionDFSData )
   {
      int newSize = maxValue(2 * newRow->memIntersectionDFSData, largestID);
      SCIP_CALL(
         SCIPreallocBlockMemoryArray(dec->env, &newRow->intersectionDFSData, newRow->memIntersectionDFSData, newSize));
      newRow->memIntersectionDFSData = newSize;
   }
   if( largestID > newRow->memColorDFSData )
   {
      int newSize = maxValue(2 * newRow->memColorDFSData, largestID);
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->colorDFSData, newRow->memColorDFSData, newSize));
      newRow->memColorDFSData = newSize;
   }
   if( largestID > newRow->memArtDFSData )
   {
      int newSize = maxValue(2 * newRow->memArtDFSData, largestID);
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->artDFSData, newRow->memArtDFSData, newSize));
      newRow->memArtDFSData = newSize;
   }

   for( int i = 0; i < newRow->memIntersectionPathDepth; ++i )
   {
      newRow->intersectionPathDepth[i] = -1;
   }

   if( largestID > newRow->memIntersectionPathDepth )
   {
      int newSize = maxValue(2 * newRow->memIntersectionPathDepth, largestID);
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->intersectionPathDepth, newRow->memIntersectionPathDepth,
                                            newSize));
      for( int i = newRow->memIntersectionPathDepth; i < newSize; ++i )
      {
         newRow->intersectionPathDepth[i] = -1;
      }
      newRow->memIntersectionPathDepth = newSize;
   }
   for( int i = 0; i < newRow->memIntersectionPathParent; ++i )
   {
      newRow->intersectionPathParent[i] = SPQR_INVALID_NODE;
   }
   if( largestID > newRow->memIntersectionPathParent )
   {
      int newSize = maxValue(2 * newRow->memIntersectionPathParent, largestID);
      SCIP_CALL(
         SCIPreallocBlockMemoryArray(dec->env, &newRow->intersectionPathParent, newRow->memIntersectionPathParent,
                                     newSize));
      for( int i = newRow->memIntersectionPathParent; i < newSize; ++i )
      {
         newRow->intersectionPathParent[i] = SPQR_INVALID_NODE;
      }
      newRow->memIntersectionPathParent = newSize;
   }

   return SCIP_OKAY;
}

static void zeroOutColors(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const spqr_node firstRemoveNode
)
{
   assert(firstRemoveNode < newRow->memNodeColors);

   newRow->nodeColors[firstRemoveNode] = UNCOLORED;
   ColorDFSCallData* data = newRow->colorDFSData;
   if( newRow->colorDFSData == NULL)
   {
      return;
   }

   data[0].node = firstRemoveNode;
   data[0].arc = getFirstNodeArc(dec, firstRemoveNode);
   if( SPQRarcIsInvalid(data[0].arc))
   {
      return;
   }

   int depth = 0;

   while( depth >= 0 )
   {
      assert(depth < newRow->memColorDFSData);
      ColorDFSCallData* callData = &data[depth];
      spqr_node head = findArcHead(dec, callData->arc);
      spqr_node tail = findArcTail(dec, callData->arc);
      spqr_node otherNode = callData->node == head ? tail : head;
      assert(otherNode < newRow->memNodeColors);
      if( newRow->nodeColors[otherNode] != UNCOLORED )
      {
         callData->arc = getNextNodeArc(dec, callData->arc, callData->node);

         newRow->nodeColors[otherNode] = UNCOLORED;
         ++depth;
         data[depth].node = otherNode;
         data[depth].arc = getFirstNodeArc(dec, otherNode);
         continue;
      }

      callData->arc = getNextNodeArc(dec, callData->arc, callData->node);
      while( depth >= 0 && data[depth].arc == getFirstNodeArc(dec, data[depth].node))
      {
         --depth;
      }
   }
}

static void cleanUpPreviousIteration(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow
)
{
   //zero out coloring information from previous check
   for( int i = 0; i < newRow->numReducedMembers; ++i )
   {
      if( SPQRnodeIsValid(newRow->reducedMembers[i].coloredNode))
      {
         zeroOutColors(dec, newRow, newRow->reducedMembers[i].coloredNode);
      }
   }

#ifndef NDEBUG
   for( int i = 0; i < newRow->memNodeColors; ++i )
   {
      assert(newRow->nodeColors[i] == UNCOLORED);
   }
#endif

   //For cut arcs: clear them from the array from previous iteration
   cut_arc_id cutArcIdx = newRow->firstOverallCutArc;
   while( cutArcIsValid(cutArcIdx))
   {
      spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
      cutArcIdx = newRow->cutArcs[cutArcIdx].nextOverall;
      newRow->isArcCut[cutArc] = FALSE;
      newRow->isArcCutReversed[cutArc] = FALSE;
   }
#ifndef NDEBUG
   for( int i = 0; i < newRow->memIsArcCut; ++i )
   {
      assert(!newRow->isArcCut[i]);
      assert(!newRow->isArcCutReversed[i]);
   }
#endif
}

static void rigidFindStarNodes(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const reduced_member_id toCheck
)
{
   //4 cases:
   //Only a single edge; both head/tail are okay => network
   //All are adjacent to a single node, but do not have it as head or tail => not network
   //All are adjacent to a single node, and have it as head or tail => network
   //Not all are adjacent to a single node => check articulation nodes
   assert(newRow->reducedMembers[toCheck].numCutArcs > 0);//calling this function otherwise is nonsensical

   cut_arc_id cutArcIdx = newRow->reducedMembers[toCheck].firstCutArc;
   spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
   spqr_node head = findArcHead(dec, cutArc);
   spqr_node tail = findArcTail(dec, cutArc);

   SCIP_Bool reverse = findArcSign(dec, cutArc).reversed != newRow->cutArcs[cutArcIdx].arcReversed;
   spqr_node cutArcsHead = reverse ? tail : head;
   spqr_node cutArcsTail = reverse ? head : tail;

   if( newRow->reducedMembers[toCheck].numCutArcs == 1 )
   {
      newRow->reducedMembers[toCheck].articulationArc = cutArc;
      newRow->reducedMembers[toCheck].splitNode = cutArcsHead;
      newRow->reducedMembers[toCheck].otherNode = cutArcsTail;
      newRow->reducedMembers[toCheck].otherIsSource = TRUE;
      newRow->reducedMembers[toCheck].allHaveCommonNode = TRUE;
      return;// Only a single cut arc
   }

   spqr_node intersectionNodes[2] = {head, tail};

   while( cutArcIsValid(newRow->cutArcs[cutArcIdx].nextMember))
   {
      cutArcIdx = newRow->cutArcs[cutArcIdx].nextMember;
      cutArc = newRow->cutArcs[cutArcIdx].arc;
      head = findArcHead(dec, cutArc);
      tail = findArcTail(dec, cutArc);
      reverse = findArcSign(dec, cutArc).reversed != newRow->cutArcs[cutArcIdx].arcReversed;
      spqr_node effectiveHead = reverse ? tail : head;
      spqr_node effectiveTail = reverse ? head : tail;
      if( effectiveHead != cutArcsHead )
      {
         cutArcsHead = SPQR_INVALID_NODE;
      }
      if( effectiveTail != cutArcsTail )
      {
         cutArcsTail = SPQR_INVALID_NODE;
      }

      //intersection between intersectionNodes and head and tail
      for( int i = 0; i < 2; ++i )
      {
         if( intersectionNodes[i] != head && intersectionNodes[i] != tail )
         {
            intersectionNodes[i] = SPQR_INVALID_NODE;
         }
      }
      if( SPQRnodeIsInvalid(intersectionNodes[0]) && SPQRnodeIsInvalid(intersectionNodes[1]))
      {
         newRow->reducedMembers[toCheck].splitNode = SPQR_INVALID_NODE;
         newRow->reducedMembers[toCheck].allHaveCommonNode = FALSE;
         return;//not all arcs are adjacent to a single node, need to check articulation nodes
      }
   }
   if( SPQRnodeIsInvalid(cutArcsHead) && SPQRnodeIsInvalid(cutArcsTail))
   {
      //All arcs adjacent to a single node, but not in same direction; not network
      newRow->remainsNetwork = FALSE;
      newRow->reducedMembers[toCheck].type = TYPE_NOT_NETWORK;
      return;
   }
   SCIP_Bool headSplittable = SPQRnodeIsValid(cutArcsHead);
   //Check if the n arcs are in a n+1 degree node; if so, the other endpoint of this non split arc is also splittable
   //By virtue of the spanning tree, this arc must be a tree arc.
   spqr_node splitNode = headSplittable ? cutArcsHead : cutArcsTail;
   newRow->reducedMembers[toCheck].splitNode = splitNode;
   newRow->reducedMembers[toCheck].otherIsSource = headSplittable;
   newRow->reducedMembers[toCheck].allHaveCommonNode = TRUE;
   if( newRow->reducedMembers[toCheck].numCutArcs == nodeDegree(dec, splitNode) - 1 )
   {
      spqr_arc firstNodeArc = getFirstNodeArc(dec, splitNode);
      spqr_arc neighbourArc = firstNodeArc;
      do
      {
         if( arcIsTree(dec, neighbourArc))
         {
            break;
         }
         neighbourArc = getNextNodeArc(dec, neighbourArc, splitNode);
      } while( neighbourArc != firstNodeArc );

      newRow->reducedMembers[toCheck].articulationArc = neighbourArc;
      spqr_arc arcHead = findArcHead(dec, neighbourArc);
      newRow->reducedMembers[toCheck].otherNode = arcHead == splitNode ? findArcTail(dec, neighbourArc) : arcHead;
   }
}

static void zeroOutColorsExceptNeighbourhood(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const spqr_node articulationNode,
   const spqr_node startRemoveNode
)
{
   COLOR_STATUS* neighbourColors = newRow->temporaryColorArray;
   assert(nodeDegree(dec, articulationNode) <= newRow->memTemporaryColorArray);
   {
      int i = 0;
      spqr_arc artFirstArc = getFirstNodeArc(dec, articulationNode);
      spqr_arc artItArc = artFirstArc;
      do
      {
         spqr_node head = findArcHead(dec, artItArc);
         spqr_node tail = findArcTail(dec, artItArc);
         spqr_node otherNode = articulationNode == head ? tail : head;
         neighbourColors[i] = newRow->nodeColors[otherNode];
         i++;
         assert(i <= nodeDegree(dec, articulationNode));
         artItArc = getNextNodeArc(dec, artItArc, articulationNode);
      } while( artItArc != artFirstArc );
   }
   zeroOutColors(dec, newRow, startRemoveNode);

   {
      int i = 0;
      spqr_arc artFirstArc = getFirstNodeArc(dec, articulationNode);
      spqr_arc artItArc = artFirstArc;
      do
      {
         spqr_node head = findArcHead(dec, artItArc);
         spqr_node tail = findArcTail(dec, artItArc);
         spqr_node otherNode = articulationNode == head ? tail : head;
         newRow->nodeColors[otherNode] = neighbourColors[i];
         i++;
         assert(i <= nodeDegree(dec, articulationNode));
         artItArc = getNextNodeArc(dec, artItArc, articulationNode);
      } while( artItArc != artFirstArc );
   }

}

//Theoretically, there is a faster algorithm, but it is quite complicated to implement.
//Thus, we stick with the 'simple' version below, which is easily fast enough in practice.
static void intersectionOfAllPaths(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const reduced_member_id toCheck,
   int* const nodeNumPaths
)
{
   int* intersectionPathDepth = newRow->intersectionPathDepth;
   spqr_node* intersectionPathParent = newRow->intersectionPathParent;

   //First do a dfs over the tree, storing all the tree-parents and depths for each node
   //TODO: maybe cache this tree and also update it so we can prevent this DFS call?

   //pick an arbitrary node as root; we just use the first cutArc here
   {
      spqr_node root = findArcHead(dec, newRow->cutArcs[newRow->reducedMembers[toCheck].firstCutArc].arc);
      DFSCallData* pathSearchCallStack = newRow->intersectionDFSData;

      assert(intersectionPathDepth[root] == -1);
      assert(intersectionPathParent[root] == SPQR_INVALID_NODE);

      int pathSearchCallStackSize = 0;

      intersectionPathDepth[root] = 0;
      intersectionPathParent[root] = SPQR_INVALID_NODE;

      pathSearchCallStack[0].node = root;
      pathSearchCallStack[0].nodeArc = getFirstNodeArc(dec, root);
      pathSearchCallStackSize++;
      while( pathSearchCallStackSize > 0 )
      {
         assert(pathSearchCallStackSize <= newRow->memIntersectionDFSData);
         DFSCallData* dfsData = &pathSearchCallStack[pathSearchCallStackSize - 1];
         //cannot be a tree arc which is its parent
         if( arcIsTree(dec, dfsData->nodeArc) &&
             ( pathSearchCallStackSize <= 1 ||
               dfsData->nodeArc != pathSearchCallStack[pathSearchCallStackSize - 2].nodeArc ))
         {
            spqr_node head = findArcHeadNoCompression(dec, dfsData->nodeArc);
            spqr_node tail = findArcTailNoCompression(dec, dfsData->nodeArc);
            spqr_node other = head == dfsData->node ? tail : head;
            assert(other != dfsData->node);

            //We go up a level: add new node to the call stack
            pathSearchCallStack[pathSearchCallStackSize].node = other;
            pathSearchCallStack[pathSearchCallStackSize].nodeArc = getFirstNodeArc(dec, other);
            //Every time a new node is discovered/added, we update its parent and depth information
            assert(intersectionPathDepth[other] == -1);
            assert(intersectionPathParent[other] == SPQR_INVALID_NODE);
            intersectionPathParent[other] = dfsData->node;
            intersectionPathDepth[other] = pathSearchCallStackSize;
            ++pathSearchCallStackSize;
            continue;
         }
         do
         {
            dfsData->nodeArc = getNextNodeArc(dec, dfsData->nodeArc, dfsData->node);
            if( dfsData->nodeArc == getFirstNodeArc(dec, dfsData->node))
            {
               --pathSearchCallStackSize;
               dfsData = &pathSearchCallStack[pathSearchCallStackSize - 1];
            } else
            {
               break;
            }
         } while( pathSearchCallStackSize > 0 );
      }
   }

   //For each cut arc, trace back both ends until they meet
   cut_arc_id cutArc = newRow->reducedMembers[toCheck].firstCutArc;
   do
   {
      spqr_arc arc = newRow->cutArcs[cutArc].arc;
      cutArc = newRow->cutArcs[cutArc].nextMember;

      //Iteratively jump up to the parents until they reach a common parent
      spqr_node source = findArcHead(dec, arc);
      spqr_node target = findArcTail(dec, arc);
      int sourceDepth = intersectionPathDepth[source];
      int targetDepth = intersectionPathDepth[target];
      nodeNumPaths[source]++;
      nodeNumPaths[target]++;

      while( sourceDepth > targetDepth )
      {
         assert(source != target);
         source = intersectionPathParent[source];
         nodeNumPaths[source]++;
         --sourceDepth;
      }
      while( targetDepth > sourceDepth )
      {
         assert(source != target);
         target = intersectionPathParent[target];
         nodeNumPaths[target]++;
         --targetDepth;
      }
      while( source != target && targetDepth >= 0 )
      {
         source = intersectionPathParent[source];
         target = intersectionPathParent[target];
         nodeNumPaths[source]++;
         nodeNumPaths[target]++;
         --targetDepth;
      }
      //In all the above, the lowest common ancestor is increased twice, so we correct for it ad-hoc
      nodeNumPaths[source]--;
      assert(SPQRnodeIsValid(source) && SPQRnodeIsValid(target));
      assert(source == target);

   } while( cutArcIsValid(cutArc));
}

static void addArticulationNode(
   SCIP_NETROWADD* newRow,
   spqr_node articulationNode
)
{
#ifndef NDEBUG
   for( int i = 0; i < newRow->numArticulationNodes; ++i )
   {
      assert(newRow->articulationNodes[i] != articulationNode);
   }
#endif
   newRow->articulationNodes[newRow->numArticulationNodes] = articulationNode;
   ++newRow->numArticulationNodes;
}

static void articulationPoints(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   ArticulationNodeInformation* nodeInfo,
   reduced_member_id reducedMember
)
{
   const SCIP_Bool* arcRemoved = newRow->isArcCut;

   int rootChildren = 0;
   spqr_node root_node = findArcHead(dec, getFirstMemberArc(dec, newRow->reducedMembers[reducedMember].member));

   ArticulationPointCallStack* callStack = newRow->artDFSData;

   int depth = 0;
   int time = 1;

   callStack[depth].arc = getFirstNodeArc(dec, root_node);
   callStack[depth].node = root_node;
   callStack[depth].parent = SPQR_INVALID_NODE;
   callStack[depth].isAP = FALSE;

   nodeInfo[root_node].low = time;
   nodeInfo[root_node].discoveryTime = time;

   while( depth >= 0 )
   {
      if( !arcRemoved[callStack[depth].arc] )
      {
         spqr_node node = callStack[depth].node;
         spqr_node head = findArcHead(dec, callStack[depth].arc);
         spqr_node tail = findArcTail(dec, callStack[depth].arc);
         spqr_node otherNode = node == head ? tail : head;
         if( otherNode != callStack[depth].parent )
         {
            if( nodeInfo[otherNode].discoveryTime == 0 )
            {
               if( depth == 0 )
               {
                  rootChildren++;
               }
               //recursive call
               ++depth;
               assert(depth < newRow->memArtDFSData);
               callStack[depth].parent = node;
               callStack[depth].node = otherNode;
               callStack[depth].arc = getFirstNodeArc(dec, otherNode);
               callStack[depth].isAP = FALSE;

               ++time;
               nodeInfo[otherNode].low = time;
               nodeInfo[otherNode].discoveryTime = time;
               continue;

            } else
            {
               nodeInfo[node].low = minValue(nodeInfo[node].low, nodeInfo[otherNode].discoveryTime);
            }
         }
      }

      while( TRUE ) /*lint !e716*/
      {
         callStack[depth].arc = getNextNodeArc(dec, callStack[depth].arc, callStack[depth].node);
         if( callStack[depth].arc != getFirstNodeArc(dec, callStack[depth].node)) break;
         --depth;
         if( depth < 0 ) break;

         spqr_node current_node = callStack[depth].node;
         spqr_node other_node = callStack[depth + 1].node;
         nodeInfo[current_node].low = minValue(nodeInfo[current_node].low,
                                          nodeInfo[other_node].low);
         if( depth != 0 &&
             !callStack[depth].isAP &&
             nodeInfo[current_node].discoveryTime <= nodeInfo[other_node].low )
         {
            addArticulationNode(newRow, current_node);
            callStack[depth].isAP = TRUE;
         }
      }
   }
   if( rootChildren > 1 )
   {
      addArticulationNode(newRow, root_node);
   }
}

static void rigidConnectedColoringRecursive(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   spqr_node articulationNode,
   spqr_node firstProcessNode,
   COLOR_STATUS firstColor,
   SCIP_Bool* isGood
)
{
   const SCIP_Bool* isArcCut = newRow->isArcCut;
   COLOR_STATUS* nodeColors = newRow->nodeColors;
   ColorDFSCallData* data = newRow->colorDFSData;

   data[0].node = firstProcessNode;
   data[0].arc = getFirstNodeArc(dec, firstProcessNode);
   newRow->nodeColors[firstProcessNode] = firstColor;

   int depth = 0;
   while( depth >= 0 )
   {
      assert(depth < newRow->memColorDFSData);
      assert(newRow->nodeColors[articulationNode] == UNCOLORED);

      ColorDFSCallData* callData = &data[depth];
      spqr_node head = findArcHead(dec, callData->arc);
      spqr_node tail = findArcTail(dec, callData->arc);
      spqr_node otherNode = callData->node == head ? tail : head;
      COLOR_STATUS currentColor = nodeColors[callData->node];
      COLOR_STATUS otherColor = nodeColors[otherNode];
      //Checks the direction of the arc; in the rest of the algorithm, we just need to check partition
      if( isArcCut[callData->arc] && currentColor != otherColor )
      {
         SCIP_Bool otherIsTail = callData->node == head;
         SCIP_Bool arcReversed = findArcSign(dec, callData->arc).reversed != newRow->isArcCutReversed[callData->arc];
         SCIP_Bool good = ( currentColor == COLOR_SOURCE ) == ( otherIsTail == arcReversed );
         if( !good )
         {
            *isGood = FALSE;
            break;
         }
      }
      if( otherNode != articulationNode )
      {
         if( otherColor == UNCOLORED )
         {
            if( isArcCut[callData->arc] )
            {
               nodeColors[otherNode] = currentColor == COLOR_SOURCE ? COLOR_SINK : COLOR_SOURCE;//reverse the colors
            } else
            {
               nodeColors[otherNode] = currentColor;
            }
            callData->arc = getNextNodeArc(dec, callData->arc, callData->node);

            depth++;
            assert(depth < newRow->memColorDFSData);
            data[depth].node = otherNode;
            data[depth].arc = getFirstNodeArc(dec, otherNode);
            continue;
         }
         if( isArcCut[callData->arc] != ( currentColor != otherColor ))
         {
            *isGood = FALSE;
            break;
         }
      }
      callData->arc = getNextNodeArc(dec, callData->arc, callData->node);
      while( depth >= 0 && data[depth].arc == getFirstNodeArc(dec, data[depth].node))
      {
         --depth;
      }
   }
}

static void rigidConnectedColoring(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const reduced_member_id reducedMember,
   const spqr_node node,
   SCIP_Bool* const isGood
)
{

   //we should only perform this function if there's more than one cut arc
   assert(newRow->reducedMembers[reducedMember].numCutArcs > 1);
#ifndef NDEBUG
   {

      spqr_member member = newRow->reducedMembers[reducedMember].member;
      spqr_arc firstArc = getFirstMemberArc(dec, member);
      spqr_arc memberArc = firstArc;
      do
      {
         assert(newRow->nodeColors[findArcHeadNoCompression(dec, memberArc)] == UNCOLORED);
         assert(newRow->nodeColors[findArcTailNoCompression(dec, memberArc)] == UNCOLORED);
         memberArc = getNextMemberArc(dec, memberArc);
      } while( firstArc != memberArc );
   }
#endif

   spqr_node firstProcessNode;
   COLOR_STATUS firstColor;
   {
      cut_arc_id cutArc = newRow->reducedMembers[reducedMember].firstCutArc;
      spqr_arc arc = newRow->cutArcs[cutArc].arc;
      assert(SPQRarcIsValid(arc));
      spqr_node head = findArcHead(dec, arc);
      spqr_node tail = findArcTail(dec, arc);
      if( findArcSign(dec, arc).reversed != newRow->cutArcs[cutArc].arcReversed )
      {
         spqr_node temp = head;
         head = tail;
         tail = temp;
      }
      if( tail != node )
      {
         firstProcessNode = tail;
         firstColor = COLOR_SOURCE;
      } else
      {
         assert(head != node);
         firstProcessNode = head;
         firstColor = COLOR_SINK;
      }
   }
   assert(SPQRnodeIsValid(firstProcessNode) && firstProcessNode != node);
   *isGood = TRUE;
   newRow->reducedMembers[reducedMember].coloredNode = firstProcessNode;
   rigidConnectedColoringRecursive(dec, newRow, node, firstProcessNode, firstColor, isGood);

   // Need to zero all colors for next attempts if we failed
   if( !( *isGood ))
   {
      zeroOutColors(dec, newRow, firstProcessNode);
      newRow->reducedMembers[reducedMember].coloredNode = SPQR_INVALID_NODE;
   } else
   {
      //Otherwise, we zero out all colors but the ones which we need
      zeroOutColorsExceptNeighbourhood(dec, newRow, node, firstProcessNode);
      newRow->reducedMembers[reducedMember].coloredNode = node;
   }
}

static spqr_node checkNeighbourColoringArticulationNode(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const spqr_node articulationNode,
   spqr_arc* const adjacentSplittingArc
)
{
   spqr_node firstSideCandidate = SPQR_INVALID_NODE;
   spqr_node secondSideCandidate = SPQR_INVALID_NODE;
   spqr_arc firstSideArc = SPQR_INVALID_ARC;
   spqr_arc secondSideArc = SPQR_INVALID_ARC;
   int numFirstSide = 0;
   int numSecondSide = 0;

   spqr_arc firstArc = getFirstNodeArc(dec, articulationNode);
   spqr_arc moveArc = firstArc;
   do
   {
      spqr_node head = findArcHead(dec, moveArc);
      spqr_node tail = findArcTail(dec, moveArc);
      spqr_node otherNode = articulationNode == head ? tail : head;
      assert(newRow->nodeColors[otherNode] != UNCOLORED);
      if(( newRow->nodeColors[otherNode] == COLOR_SOURCE ) != newRow->isArcCut[moveArc] )
      {
         if( numFirstSide == 0 && arcIsTree(dec, moveArc))
         {
            firstSideCandidate = otherNode;
            firstSideArc = moveArc;
         } else if( numFirstSide > 0 )
         {
            firstSideCandidate = SPQR_INVALID_NODE;
         }
         ++numFirstSide;
      } else
      {
         if( numSecondSide == 0 && arcIsTree(dec, moveArc))
         {
            secondSideCandidate = otherNode;
            secondSideArc = moveArc;
         } else if( numSecondSide > 0 )
         {
            secondSideCandidate = SPQR_INVALID_NODE;
         }
         ++numSecondSide;
      }
      moveArc = getNextNodeArc(dec, moveArc, articulationNode);
   } while( moveArc != firstArc );

   if( numFirstSide == 1 )
   {
      *adjacentSplittingArc = firstSideArc;
      return firstSideCandidate;
   } else if( numSecondSide == 1 )
   {
      *adjacentSplittingArc = secondSideArc;
      return secondSideCandidate;
   }
   return SPQR_INVALID_NODE;
}


static void rigidGetSplittableArticulationPointsOnPath(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const reduced_member_id toCheck,
   spqr_node firstNode,
   spqr_node secondNode
)
{
   assert(newRow->reducedMembers[toCheck].numCutArcs > 1);

   int totalNumNodes = getNumNodes(dec);
   int* nodeNumPaths = newRow->crossingPathCount;

   for( int i = 0; i < totalNumNodes; ++i )
   {
      nodeNumPaths[i] = 0;
   }

   intersectionOfAllPaths(dec, newRow, toCheck, nodeNumPaths);

   newRow->numArticulationNodes = 0;

   ArticulationNodeInformation* artNodeInfo = newRow->articulationNodeSearchInfo;
   //TODO: ugly; we clean up over all decomposition nodes for every component
   //clean up can not easily be done in the search, unfortunately
   for( int i = 0; i < totalNumNodes; ++i )
   {
      artNodeInfo[i].low = 0;
      artNodeInfo[i].discoveryTime = 0;
   }

   articulationPoints(dec, newRow, artNodeInfo, toCheck);

   int numCutArcs = newRow->reducedMembers[toCheck].numCutArcs;
   for( int i = 0; i < newRow->numArticulationNodes; i++ )
   {
      spqr_node articulationNode = newRow->articulationNodes[i];
      assert(nodeIsRepresentative(dec, articulationNode));
      SCIP_Bool isOnPath = nodeNumPaths[articulationNode] == numCutArcs;
      if( isOnPath &&
          (( SPQRnodeIsInvalid(firstNode) && SPQRnodeIsInvalid(secondNode)) || articulationNode == firstNode ||
           articulationNode == secondNode ))
      {
         SCIP_Bool isGood = TRUE;
         rigidConnectedColoring(dec, newRow, toCheck, articulationNode, &isGood);
         if( !isGood )
         {
            continue;
         }
         newRow->reducedMembers[toCheck].splitNode = articulationNode;

         //Once we have found one node, we can (possibly) find another by looking at the coloring of the neighbourhood of it.

         spqr_arc adjacentSplittingArc = SPQR_INVALID_ARC;
         spqr_node adjacentSplittingNode = checkNeighbourColoringArticulationNode(dec, newRow, articulationNode,
                                                                                  &adjacentSplittingArc);

         //If the neighbour-coloring works, we still need to check that the adjacent node
         //is also an articulation node
         if( SPQRnodeIsValid(adjacentSplittingNode) &&
             (( SPQRnodeIsInvalid(firstNode) && SPQRnodeIsInvalid(secondNode)) || adjacentSplittingNode == firstNode ||
              adjacentSplittingNode == secondNode ))
         {
            SCIP_Bool isArticulationNode = FALSE;
            for( int j = 0; j < newRow->numArticulationNodes; ++j )
            {
               if( newRow->articulationNodes[j] == adjacentSplittingNode )
               {
                  isArticulationNode = TRUE;
                  break;
               }
            }
            if( isArticulationNode )
            {
               newRow->reducedMembers[toCheck].articulationArc = adjacentSplittingArc;
               newRow->reducedMembers[toCheck].otherNode = adjacentSplittingNode;
               newRow->reducedMembers[toCheck].otherIsSource =
                  newRow->nodeColors[adjacentSplittingNode] == COLOR_SOURCE;

               //Cleaning up the colors
               {
                  spqr_arc firstNodeArc = getFirstNodeArc(dec, articulationNode);
                  spqr_arc itArc = firstNodeArc;
                  do
                  {
                     spqr_node head = findArcHead(dec, itArc);
                     spqr_node tail = findArcTail(dec, itArc);
                     spqr_node otherNode = articulationNode == head ? tail : head;
                     newRow->nodeColors[otherNode] = UNCOLORED;
                     itArc = getNextNodeArc(dec, itArc, articulationNode);
                  } while( itArc != firstNodeArc );
               }
            }
         }
         return;
      }
   }

   newRow->reducedMembers[toCheck].type = TYPE_NOT_NETWORK;
   newRow->remainsNetwork = FALSE;
}

static void determineParallelType(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const reduced_member_id toCheckMember,
   const spqr_arc markerToOther,
   const reduced_member_id otherMember,
   const spqr_arc markerToCheck
)
{

   SPQRRowReducedMember* member = &newRow->reducedMembers[toCheckMember];
   assert(cutArcIsValid(member->firstCutArc));

   SCIP_Bool good = TRUE;
   SCIP_Bool isReversed = TRUE;
   int countedCutArcs = 0;
   for( cut_arc_id cutArc = member->firstCutArc; cutArcIsValid(cutArc);
        cutArc = newRow->cutArcs[cutArc].nextMember )
   {
      spqr_arc arc = newRow->cutArcs[cutArc].arc;
      SCIP_Bool arcIsReversed = arcIsReversedNonRigid(dec, arc) != newRow->cutArcs[cutArc].arcReversed;
      if( countedCutArcs == 0 )
      {
         isReversed = arcIsReversed;
      } else if( arcIsReversed != isReversed )
      {
         good = FALSE;
         break;
      }
      ++countedCutArcs;
   }
   if( !good )
   {
      member->type = TYPE_NOT_NETWORK;
      newRow->remainsNetwork = FALSE;
      return;
   }
   //we can only propagate if the marker arc is a tree arc and all other arcs are cut
   if( !arcIsTree(dec, markerToOther) ||
       countedCutArcs != ( getNumMemberArcs(dec, newRow->reducedMembers[toCheckMember].member) - 1 ))
   {
      //In all other cases, the bond can be split so that the result will be okay!
      newRow->reducedMembers[toCheckMember].type = TYPE_MERGED;
   } else
   {
      SCIP_Bool markerIsReversed = arcIsReversedNonRigid(dec, markerToOther);
      createCutArc(dec, newRow, markerToCheck, otherMember, markerIsReversed != isReversed);
      newRow->reducedMembers[toCheckMember].type = TYPE_PROPAGATED;
   }
}

static void determineSeriesType(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const reduced_member_id toCheckMember,
   const spqr_arc markerToOther,
   const reduced_member_id otherMember,
   const spqr_arc markerToCheck
)
{
   //Propagation only calls this function if the arc is tree already, so we do not check it here.
   assert(arcIsTree(dec, markerToOther));
   assert(newRow->reducedMembers[toCheckMember].numCutArcs == 1);
   assert(cutArcIsValid(newRow->reducedMembers[toCheckMember].firstCutArc));
   spqr_arc cutArc = newRow->reducedMembers[toCheckMember].firstCutArc;
   spqr_arc arc = newRow->cutArcs[cutArc].arc;
   newRow->reducedMembers[toCheckMember].type = TYPE_PROPAGATED;
   createCutArc(dec, newRow, markerToCheck, otherMember,
                ( arcIsReversedNonRigid(dec, arc) == arcIsReversedNonRigid(dec, markerToOther)) !=
                newRow->cutArcs[cutArc].arcReversed);
}

static void determineRigidType(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const reduced_member_id toCheckMember,
   const spqr_arc markerToOther,
   const reduced_member_id otherMember,
   const spqr_arc markerToCheck
)
{
   assert(newRow->reducedMembers[toCheckMember].numCutArcs >
          0);//Checking for propagation only makes sense if there is at least one cut arc

   rigidFindStarNodes(dec, newRow, toCheckMember);
   if( !newRow->remainsNetwork )
   {
      return;
   }
   spqr_node markerHead = findArcHead(dec, markerToOther);
   spqr_node markerTail = findArcTail(dec, markerToOther);
   if( findArcSign(dec, markerToOther).reversed )
   {
      spqr_node temp = markerHead;
      markerHead = markerTail;
      markerTail = temp;
   }
   if( SPQRnodeIsInvalid(newRow->reducedMembers[toCheckMember].splitNode))
   {
      //not a star => attempt to find splittable nodes using articulation node algorithms
      //We save some work by telling the methods that only the marker nodes should be checked
      rigidGetSplittableArticulationPointsOnPath(dec, newRow, toCheckMember, markerHead, markerTail);
   }
   if( !newRow->remainsNetwork )
   {
      return;
   }


   if( SPQRarcIsValid(newRow->reducedMembers[toCheckMember].articulationArc) &&
       newRow->reducedMembers[toCheckMember].articulationArc == markerToOther )
   {
      newRow->reducedMembers[toCheckMember].type = TYPE_PROPAGATED;
      SCIP_Bool reverse = ( markerHead == newRow->reducedMembers[toCheckMember].splitNode ) !=
                          newRow->reducedMembers[toCheckMember].otherIsSource;

      createCutArc(dec, newRow, markerToCheck, otherMember, reverse);
   } else if( newRow->reducedMembers[toCheckMember].splitNode == markerHead ||
              newRow->reducedMembers[toCheckMember].splitNode == markerTail )
   {
      newRow->reducedMembers[toCheckMember].type = TYPE_MERGED;
   } else if( SPQRarcIsValid(newRow->reducedMembers[toCheckMember].articulationArc) &&
              ( newRow->reducedMembers[toCheckMember].otherNode == markerHead ||
                newRow->reducedMembers[toCheckMember].otherNode == markerTail ))
   {
      newRow->reducedMembers[toCheckMember].type = TYPE_MERGED;
   } else
   {
      //Found source or sinks, but not adjacent to the marker
      newRow->reducedMembers[toCheckMember].type = TYPE_NOT_NETWORK;
      newRow->remainsNetwork = FALSE;
   }
}

static RowReducedMemberType determineType(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const reduced_member_id toCheckMember,
   const spqr_arc markerToOther,
   const reduced_member_id otherMember,
   const spqr_arc markerToCheck
)
{
   assert(newRow->reducedMembers[toCheckMember].type == TYPE_UNDETERMINED);
   switch( getMemberType(dec, newRow->reducedMembers[toCheckMember].member))
   {
      case SPQR_MEMBERTYPE_RIGID:
      {
         determineRigidType(dec, newRow, toCheckMember, markerToOther, otherMember, markerToCheck);
         break;
      }
      case SPQR_MEMBERTYPE_LOOP:
      case SPQR_MEMBERTYPE_PARALLEL:
      {
         determineParallelType(dec, newRow, toCheckMember, markerToOther, otherMember, markerToCheck);
         break;
      }
      case SPQR_MEMBERTYPE_SERIES:
      {
         determineSeriesType(dec, newRow, toCheckMember, markerToOther, otherMember, markerToCheck);
         break;
      }
      case SPQR_MEMBERTYPE_UNASSIGNED:
      default:
         newRow->reducedMembers[toCheckMember].type = TYPE_NOT_NETWORK;
         assert(FALSE);
   }

   return newRow->reducedMembers[toCheckMember].type;
}

static void propagateComponents(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow
)
{
   int leafArrayIndex = 0;

   reduced_member_id leaf;
   reduced_member_id next;

   while( leafArrayIndex != newRow->numLeafMembers )
   {
      leaf = newRow->leafMembers[leafArrayIndex];
      //next is invalid if the member is not in the reduced decomposition.
      next = newRow->reducedMembers[leaf].parent;
      spqr_arc parentMarker = markerToParent(dec, newRow->reducedMembers[leaf].member);
      if( next != INVALID_REDUCED_MEMBER && arcIsTree(dec, parentMarker))
      {
         assert(reducedMemberIsValid(next));
         assert(SPQRarcIsValid(parentMarker));
         RowReducedMemberType type = determineType(dec, newRow, leaf, parentMarker, next,
                                                   markerOfParent(dec, newRow->reducedMembers[leaf].member));
         if( type == TYPE_PROPAGATED )
         {
            ++newRow->reducedMembers[next].numPropagatedChildren;
            if( newRow->reducedMembers[next].numPropagatedChildren == newRow->reducedMembers[next].numChildren )
            {
               newRow->leafMembers[leafArrayIndex] = next;
            } else
            {
               ++leafArrayIndex;
            }
         } else if( type == TYPE_NOT_NETWORK )
         {
            return;
         } else
         {
            assert(type == TYPE_MERGED);
            ++leafArrayIndex;
         }
      } else
      {
         ++leafArrayIndex;
      }
   }

   for( int j = 0; j < newRow->numReducedComponents; ++j )
   {
      //The reduced root might be a leaf as well: we propagate it last
      reduced_member_id root = newRow->reducedComponents[j].root;

      while( TRUE ) /*lint !e716*/
      {
         if( newRow->reducedMembers[root].numPropagatedChildren == newRow->reducedMembers[root].numChildren - 1 )
         {
            //TODO: bit ugly, have to do a linear search for the child
            reduced_member_id child = INVALID_REDUCED_MEMBER;
            spqr_arc markerToChild = SPQR_INVALID_ARC;
            for( children_idx i = newRow->reducedMembers[root].firstChild;
                 i < newRow->reducedMembers[root].firstChild + newRow->reducedMembers[root].numChildren; ++i )
            {
               reduced_member_id childReduced = newRow->childrenStorage[i];
               if( newRow->reducedMembers[childReduced].type != TYPE_PROPAGATED )
               {
                  child = childReduced;
                  markerToChild = markerOfParent(dec, newRow->reducedMembers[child].member);
                  break;
               }
            }
            assert(SPQRmemberIsValid(newRow->reducedMembers[child].member));
            assert(SPQRarcIsValid(markerToChild));
            if( !arcIsTree(dec, markerToChild))
            {
               break;
            }
            RowReducedMemberType type = determineType(dec, newRow, root, markerToChild, child,
                                                      markerToParent(dec, newRow->reducedMembers[child].member));
            if( type == TYPE_PROPAGATED )
            {
               root = child;
            } else if( type == TYPE_NOT_NETWORK )
            {
               return;
            } else
            {
               break;
            }
         } else
         {
            break;
         }
      }
      newRow->reducedComponents[j].root = root;
      newRow->reducedMembers[root].parent = INVALID_REDUCED_MEMBER;
   }
}

typedef struct
{
   spqr_node first;
   spqr_node second;
} Nodes;

static Nodes
rigidDetermineCandidateNodesFromAdjacentComponents(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const reduced_member_id toCheck
)
{

   Nodes pair;
   pair.first = SPQR_INVALID_NODE;
   pair.second = SPQR_INVALID_NODE;

   //take union of children's arcs nodes to find one or two candidates
   for( children_idx i = newRow->reducedMembers[toCheck].firstChild;
        i < newRow->reducedMembers[toCheck].firstChild + newRow->reducedMembers[toCheck].numChildren; ++i )
   {
      reduced_member_id reducedChild = newRow->childrenStorage[i];
      if( newRow->reducedMembers[reducedChild].type != TYPE_PROPAGATED )
      {
         spqr_arc arc = markerOfParent(dec, newRow->reducedMembers[reducedChild].member);
         spqr_node head = findArcHead(dec, arc);
         spqr_node tail = findArcTail(dec, arc);
         if( SPQRnodeIsInvalid(pair.first) && SPQRnodeIsInvalid(pair.second))
         {
            pair.first = head;
            pair.second = tail;
         } else
         {
            if( pair.first != head && pair.first != tail )
            {
               pair.first = SPQR_INVALID_NODE;
            }
            if( pair.second != head && pair.second != tail )
            {
               pair.second = SPQR_INVALID_NODE;
            }
         }
         if( SPQRnodeIsInvalid(pair.first) && SPQRnodeIsInvalid(pair.second))
         {
            return pair;
         }
      }
   }
   if( reducedMemberIsValid(newRow->reducedMembers[toCheck].parent) &&
       newRow->reducedMembers[newRow->reducedMembers[toCheck].parent].type != TYPE_PROPAGATED )
   {

      spqr_arc arc = markerToParent(dec, newRow->reducedMembers[toCheck].member);
      spqr_node head = findArcHead(dec, arc);
      spqr_node tail = findArcTail(dec, arc);
      if( SPQRnodeIsInvalid(pair.first) && SPQRnodeIsInvalid(pair.second))
      {
         pair.first = head;
         pair.second = tail;
      } else
      {
         if( pair.first != head && pair.first != tail )
         {
            pair.first = SPQR_INVALID_NODE;
         }
         if( pair.second != head && pair.second != tail )
         {
            pair.second = SPQR_INVALID_NODE;
         }
      }
   }
   return pair;
}

static SCIP_RETCODE allocateTreeSearchMemory(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow
)
{
   int necessarySpace = newRow->numReducedMembers;
   if( necessarySpace > newRow->memMergeTreeCallData )
   {
      int updatedSize = maxValue(2 * newRow->memMergeTreeCallData, necessarySpace);
      SCIP_CALL(SCIPreallocBlockMemoryArray(dec->env, &newRow->mergeTreeCallData, newRow->memMergeTreeCallData,
                                            updatedSize));
      newRow->memMergeTreeCallData = updatedSize;
   }
   return SCIP_OKAY;
}

static void determineSingleRowRigidType(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id reducedMember
)
{
   assert(newRow->reducedMembers[reducedMember].numCutArcs >
          0);//Checking for propagation only makes sense if there is at least one cut arc
   rigidFindStarNodes(dec, newRow, reducedMember);
   if( !newRow->remainsNetwork )
   {
      return;
   }

   if( SPQRnodeIsInvalid(newRow->reducedMembers[reducedMember].splitNode))
   {
      //not a star => attempt to find splittable nodes using articulation node algorithms
      rigidGetSplittableArticulationPointsOnPath(dec, newRow, reducedMember, SPQR_INVALID_NODE, SPQR_INVALID_NODE);
   }
   if( SPQRnodeIsValid(newRow->reducedMembers[reducedMember].splitNode))
   {
      newRow->reducedMembers[reducedMember].type = TYPE_MERGED;
   } else
   {
      newRow->reducedMembers[reducedMember].type = TYPE_NOT_NETWORK;
      newRow->remainsNetwork = FALSE;
   }
}

static void determineSingleParallelType(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id reducedMember
   )
{
   SPQRRowReducedMember* redMember = &newRow->reducedMembers[reducedMember];
   assert(cutArcIsValid(redMember->firstCutArc));

   SCIP_Bool good = TRUE;
   SCIP_Bool isReversed = TRUE;
   int countedCutArcs = 0;
   for( cut_arc_id cutArc = redMember->firstCutArc; cutArcIsValid(cutArc);
        cutArc = newRow->cutArcs[cutArc].nextMember )
   {
      spqr_arc arc = newRow->cutArcs[cutArc].arc;
      SCIP_Bool arcIsReversed = arcIsReversedNonRigid(dec, arc) != newRow->cutArcs[cutArc].arcReversed;
      if( countedCutArcs == 0 )
      {
         isReversed = arcIsReversed;
      } else if( arcIsReversed != isReversed )
      {
         good = FALSE;
         break;
      }
      ++countedCutArcs;
   }
   if( !good )
   {
      redMember->type = TYPE_NOT_NETWORK;
      newRow->remainsNetwork = FALSE;
   } else
   {
      redMember->type = TYPE_MERGED;
   }
}

static void determineSingleSeriesType(
   SCIP_NETROWADD* newRow,
   reduced_member_id reducedMember
)
{
   assert(newRow->reducedMembers[reducedMember].numCutArcs == 1);
   assert(cutArcIsValid(newRow->reducedMembers[reducedMember].firstCutArc));
   newRow->reducedMembers[reducedMember].type = TYPE_MERGED;
}

static spqr_node determineAndColorSplitNode(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id id,
   spqr_node candidateOne,
   spqr_node candidateTwo
)
{
   if( SPQRnodeIsInvalid(newRow->reducedMembers[id].splitNode))
   {
      return SPQR_INVALID_NODE;
   }
   spqr_node splitNode = newRow->reducedMembers[id].splitNode;
   if( splitNode == candidateOne || splitNode == candidateTwo )
   {
      if( newRow->reducedMembers[id].allHaveCommonNode )
      {
         spqr_arc firstNodeArc = getFirstNodeArc(dec, splitNode);
         spqr_arc iterArc = firstNodeArc;
         COLOR_STATUS color = newRow->reducedMembers[id].otherIsSource ? COLOR_SOURCE : COLOR_SINK;
         do
         {
            spqr_node head = findArcHead(dec, iterArc);
            spqr_node other = head == splitNode ? findArcTail(dec, iterArc) : head;
            newRow->nodeColors[other] = color;
            iterArc = getNextNodeArc(dec, iterArc, splitNode);
         } while( iterArc != firstNodeArc );
         newRow->reducedMembers[id].coloredNode = splitNode;
      }
      return splitNode;
   }
   splitNode = newRow->reducedMembers[id].otherNode;
   if( SPQRnodeIsInvalid(splitNode) || ( splitNode != candidateOne && splitNode != candidateTwo ))
   {
      return SPQR_INVALID_NODE;
   }
   assert(SPQRarcIsValid(newRow->reducedMembers[id].articulationArc));
   if( newRow->reducedMembers[id].allHaveCommonNode )
   {
      spqr_arc firstNodeArc = getFirstNodeArc(dec, splitNode);
      spqr_arc iterArc = firstNodeArc;
      COLOR_STATUS color;
      if( newRow->reducedMembers[id].numCutArcs == 1 )
      {
         color = newRow->reducedMembers[id].otherIsSource ? COLOR_SINK : COLOR_SOURCE;
      } else
      {
         color = newRow->reducedMembers[id].otherIsSource ? COLOR_SOURCE : COLOR_SINK;
      }
      do
      {
         spqr_node head = findArcHead(dec, iterArc);
         spqr_node other = head == splitNode ? findArcTail(dec, iterArc) : head;
         newRow->nodeColors[other] = color;
         iterArc = getNextNodeArc(dec, iterArc, splitNode);
      } while( iterArc != firstNodeArc );
      newRow->nodeColors[newRow->reducedMembers[id].splitNode] = newRow->reducedMembers[id].otherIsSource ? COLOR_SINK
                                                                                                          : COLOR_SOURCE;

   } else
   {
      COLOR_STATUS splitColor = newRow->nodeColors[splitNode];

      spqr_arc firstNodeArc = getFirstNodeArc(dec, splitNode);
      spqr_arc iterArc = firstNodeArc;
      do
      {
         spqr_node head = findArcHead(dec, iterArc);
         spqr_node other = head == splitNode ? findArcTail(dec, iterArc) : head;
         newRow->nodeColors[other] = splitColor;
         iterArc = getNextNodeArc(dec, iterArc, splitNode);
      } while( iterArc != firstNodeArc );
      newRow->nodeColors[newRow->reducedMembers[id].splitNode] = splitColor == COLOR_SOURCE ? COLOR_SINK : COLOR_SOURCE;
      newRow->nodeColors[splitNode] = UNCOLORED;
   }

   newRow->reducedMembers[id].coloredNode = splitNode;
   return splitNode;
}

static void determineSplitTypeFirstLeaf(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id reducedId
)
{
   spqr_member member = newRow->reducedMembers[reducedId].member;
   SPQRMemberType type = getMemberType(dec, member);
   assert(type == SPQR_MEMBERTYPE_PARALLEL || type == SPQR_MEMBERTYPE_RIGID);
   assert(cutArcIsValid(newRow->reducedMembers[reducedId].firstCutArc));
   SPQRRowReducedMember* redMember = &newRow->reducedMembers[reducedId];

   if( type == SPQR_MEMBERTYPE_PARALLEL )
   {
      //TODO: duplicate-ish

      SCIP_Bool good = TRUE;
      SCIP_Bool isReversed = TRUE;
      int countedCutArcs = 0;
      for( cut_arc_id cutArc = redMember->firstCutArc; cutArcIsValid(cutArc);
           cutArc = newRow->cutArcs[cutArc].nextMember )
      {
         spqr_arc arc = newRow->cutArcs[cutArc].arc;
         SCIP_Bool arcIsReversed = arcIsReversedNonRigid(dec, arc) != newRow->cutArcs[cutArc].arcReversed;
         if( countedCutArcs == 0 )
         {
            isReversed = arcIsReversed;
         } else if( arcIsReversed != isReversed )
         {
            good = FALSE;
            break;
         }
         ++countedCutArcs;
      }
      if( !good )
      {
         redMember->type = TYPE_NOT_NETWORK;
         newRow->remainsNetwork = FALSE;
      } else
      {
         spqr_arc marker = markerToParent(dec, member);
         redMember->type = TYPE_MERGED;
         redMember->splitArc = marker;
         redMember->splitHead = TRUE;
         redMember->otherIsSource = arcIsReversedNonRigid(dec, marker) == isReversed;
      }
      return;
   }
   assert(type == SPQR_MEMBERTYPE_RIGID);

   spqr_arc marker = markerToParent(dec, member);
   spqr_node markerHead = findArcHead(dec, marker);
   spqr_node markerTail = findArcTail(dec, marker);
   if( findArcSign(dec, marker).reversed )
   {
      spqr_node temp = markerHead;
      markerHead = markerTail;
      markerTail = temp;
   }

   if( !SPQRnodeIsValid(newRow->reducedMembers[reducedId].splitNode))
   {
      assert(newRow->reducedMembers[reducedId].numCutArcs >
             0);//Checking for propagation only makes sense if there is at least one cut arc

      rigidFindStarNodes(dec, newRow, reducedId);
      if( !newRow->remainsNetwork )
      {
         return;
      }

      if( SPQRnodeIsInvalid(newRow->reducedMembers[reducedId].splitNode))
      {
         //not a star => attempt to find splittable nodes using articulation node algorithms
         //We save some work by telling the methods that only the marker nodes should be checked
         rigidGetSplittableArticulationPointsOnPath(dec, newRow, reducedId, markerHead, markerTail);
      }
      if( !newRow->remainsNetwork )
      {
         return;
      }
      if( SPQRnodeIsInvalid(newRow->reducedMembers[reducedId].splitNode))
      {
         redMember->type = TYPE_NOT_NETWORK;
         newRow->remainsNetwork = FALSE;
         return;
      }
   }

   spqr_node splitNode = newRow->reducedMembers[reducedId].splitNode;
#ifndef NDEBUG
   spqr_node otherNode = newRow->reducedMembers[reducedId].otherNode;
#endif
   //We cannot have both splittable (should have been propagated)
   assert(!(( splitNode == markerTail || splitNode == markerHead ) &&
            ( otherNode == markerTail || otherNode == markerHead )));

   splitNode = determineAndColorSplitNode(dec, newRow, reducedId, markerHead, markerTail);
   if( SPQRnodeIsInvalid(splitNode))
   {
      redMember->type = TYPE_NOT_NETWORK;
      newRow->remainsNetwork = FALSE;
      return;
   }
   assert(splitNode == markerHead || splitNode == markerTail);

   newRow->reducedMembers[reducedId].splitNode = splitNode;
   newRow->reducedMembers[reducedId].willBeReversed = FALSE;
   redMember->type = TYPE_MERGED;
}

typedef struct
{
   SCIP_Bool headSplit;
   SCIP_Bool otherIsSource;
} SplitOrientation;

static SplitOrientation getRelativeOrientationRigid(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id reducedId,
   spqr_arc arcToNext
)
{
   assert(findArcMemberNoCompression(dec, arcToNext) == newRow->reducedMembers[reducedId].member);
   assert(SPQRnodeIsValid(newRow->reducedMembers[reducedId].splitNode));

   SplitOrientation orientation;
   if( newRow->reducedMembers[reducedId].numCutArcs == 0 )
   {
      spqr_node splitNode = newRow->reducedMembers[reducedId].splitNode;
      spqr_node head = findEffectiveArcHead(dec, arcToNext);

      assert(head == splitNode || splitNode == findEffectiveArcTailNoCompression(dec, arcToNext));

      orientation.headSplit = newRow->reducedMembers[reducedId].willBeReversed == ( head != splitNode );
      orientation.otherIsSource = newRow->reducedMembers[reducedId].otherIsSource;
      return orientation;
   }
   spqr_node splitNode = newRow->reducedMembers[reducedId].splitNode;
   spqr_node arcHead = findArcHead(dec, arcToNext);
   spqr_node arcTail = findArcTail(dec, arcToNext);
   if( findArcSign(dec, arcToNext).reversed )
   {
      spqr_node temp = arcHead;
      arcHead = arcTail;
      arcTail = temp;
   }
   assert(arcHead == splitNode || arcTail == splitNode);
   spqr_node other = arcHead == splitNode ? arcTail : arcHead;

   assert(newRow->nodeColors[other] == COLOR_SOURCE || newRow->nodeColors[other] == COLOR_SINK);

   if( newRow->reducedMembers[reducedId].willBeReversed )
   {
      orientation.headSplit = arcHead != splitNode;
      orientation.otherIsSource = newRow->nodeColors[other] != COLOR_SOURCE;
   } else
   {
      orientation.headSplit = arcHead == splitNode;
      orientation.otherIsSource = newRow->nodeColors[other] == COLOR_SOURCE;
   }
   return orientation;
}

static SplitOrientation getRelativeOrientationParallel(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id reducedId,
   spqr_arc arcToNext
)
{
   assert(findArcMemberNoCompression(dec, arcToNext) == newRow->reducedMembers[reducedId].member);
   assert(SPQRarcIsValid(newRow->reducedMembers[reducedId].splitArc) && SPQRarcIsValid(arcToNext));
   SplitOrientation orientation;
   orientation.otherIsSource = newRow->reducedMembers[reducedId].otherIsSource;
   if( arcIsReversedNonRigid(dec, arcToNext) == arcIsReversedNonRigid(dec, newRow->reducedMembers[reducedId].splitArc))
   {
      orientation.headSplit = newRow->reducedMembers[reducedId].splitHead;
   } else
   {
      orientation.headSplit = !newRow->reducedMembers[reducedId].splitHead;
   }
   return orientation;
}

static SplitOrientation getRelativeOrientationSeries(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id reducedId,
   spqr_arc arcToNext
)
{
   assert(findArcMemberNoCompression(dec, arcToNext) == newRow->reducedMembers[reducedId].member);
   assert(SPQRarcIsValid(newRow->reducedMembers[reducedId].splitArc) && SPQRarcIsValid(arcToNext));
   SplitOrientation orientation;

   orientation.otherIsSource = newRow->reducedMembers[reducedId].otherIsSource;
   if( arcIsReversedNonRigid(dec, arcToNext) == arcIsReversedNonRigid(dec, newRow->reducedMembers[reducedId].splitArc))
   {
      orientation.headSplit = !newRow->reducedMembers[reducedId].splitHead;
   } else
   {
      orientation.headSplit = newRow->reducedMembers[reducedId].splitHead;
   }
   return orientation;
}

static SplitOrientation getRelativeOrientation(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id reducedId,
   spqr_arc arcToNext
)
{
   switch( getMemberType(dec, newRow->reducedMembers[reducedId].member) )
   {
      case SPQR_MEMBERTYPE_RIGID:
         return getRelativeOrientationRigid(dec, newRow, reducedId, arcToNext);
      case SPQR_MEMBERTYPE_PARALLEL:
         return getRelativeOrientationParallel(dec, newRow, reducedId, arcToNext);
      case SPQR_MEMBERTYPE_SERIES:
         return getRelativeOrientationSeries(dec, newRow, reducedId, arcToNext);
      case SPQR_MEMBERTYPE_LOOP:
      case SPQR_MEMBERTYPE_UNASSIGNED:
      default:
      {
         assert(FALSE);
         SplitOrientation orientation; /*lint !e527*/
         orientation.headSplit = FALSE; /*lint !e527*/
         orientation.otherIsSource = FALSE;
         return orientation;
      }
   }
}

static void determineSplitTypeSeries(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id reducedId,
   spqr_arc marker,
   SplitOrientation previousOrientation
)
{
   int numAdjacentMembers =
      newRow->reducedMembers[reducedId].numChildren - newRow->reducedMembers[reducedId].numPropagatedChildren;
   if( reducedMemberIsValid(newRow->reducedMembers[reducedId].parent) &&
       newRow->reducedMembers[newRow->reducedMembers[reducedId].parent].type != TYPE_PROPAGATED )
   {
      ++numAdjacentMembers;
   }
   assert(numAdjacentMembers > 1);
   if( numAdjacentMembers > 2 )
   {
      newRow->remainsNetwork = FALSE;
      newRow->reducedMembers[reducedId].type = TYPE_NOT_NETWORK;
      return;
   }
   cut_arc_id cutArc = newRow->reducedMembers[reducedId].firstCutArc;
   if( cutArcIsValid(cutArc))
   {
      spqr_arc arc = newRow->cutArcs[cutArc].arc;
      SCIP_Bool good = ((( arcIsReversedNonRigid(dec, arc) == arcIsReversedNonRigid(dec, marker)) ==
                         newRow->cutArcs[cutArc].arcReversed ) == previousOrientation.headSplit ) ==
                       previousOrientation.otherIsSource;
      if( !good )
      {
         newRow->remainsNetwork = FALSE;
         newRow->reducedMembers[reducedId].type = TYPE_NOT_NETWORK;
         return;
      }
      newRow->reducedMembers[reducedId].splitArc = marker;
      newRow->reducedMembers[reducedId].splitHead = previousOrientation.headSplit;
      newRow->reducedMembers[reducedId].otherIsSource = !previousOrientation.otherIsSource;
      newRow->reducedMembers[reducedId].type = TYPE_MERGED;
      return;
   }

   newRow->reducedMembers[reducedId].splitArc = marker;
   newRow->reducedMembers[reducedId].splitHead = previousOrientation.headSplit;
   newRow->reducedMembers[reducedId].otherIsSource = previousOrientation.otherIsSource;
   newRow->reducedMembers[reducedId].type = TYPE_MERGED;
}

static void determineSplitTypeParallel(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id reducedId,
   spqr_arc marker,
   SplitOrientation previousOrientation
)
{
   SPQRRowReducedMember* redMember = &newRow->reducedMembers[reducedId];

   SCIP_Bool good = TRUE;
   SCIP_Bool isReversed = TRUE;
   int countedCutArcs = 0;
   for( cut_arc_id cutArc = redMember->firstCutArc; cutArcIsValid(cutArc);
        cutArc = newRow->cutArcs[cutArc].nextMember )
   {
      spqr_arc arc = newRow->cutArcs[cutArc].arc;
      SCIP_Bool arcIsReversed = arcIsReversedNonRigid(dec, arc) != newRow->cutArcs[cutArc].arcReversed;
      if( countedCutArcs == 0 )
      {
         isReversed = arcIsReversed;
      } else if( arcIsReversed != isReversed )
      {
         good = FALSE;
         break;
      }
      ++countedCutArcs;
   }
   if( !good )
   {
      redMember->type = TYPE_NOT_NETWORK;
      newRow->remainsNetwork = FALSE;
      return;
   }
   if( countedCutArcs == 0 )
   {
      redMember->splitArc = marker;
      redMember->splitHead = previousOrientation.headSplit;
      redMember->otherIsSource = previousOrientation.otherIsSource;
      redMember->type = TYPE_MERGED;
      return;
   }
   SCIP_Bool isHeadSourceOrTailTarget = previousOrientation.headSplit == previousOrientation.otherIsSource;
   SCIP_Bool isOkay = isHeadSourceOrTailTarget == ( isReversed == arcIsReversedNonRigid(dec, marker));
   if( !isOkay )
   {
      redMember->type = TYPE_NOT_NETWORK;
      newRow->remainsNetwork = FALSE;
      return;
   }
   redMember->splitArc = marker;
   redMember->splitHead = previousOrientation.headSplit;
   redMember->otherIsSource = previousOrientation.otherIsSource;
   redMember->type = TYPE_MERGED;
}

static void determineSplitTypeRigid(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id reducedId,
   spqr_member member,
   spqr_arc marker,
   SplitOrientation previousOrientation
)
{
   assert(dec);
   assert(newRow);
   assert(getMemberType(dec, member) == SPQR_MEMBERTYPE_RIGID);

   Nodes nodes = rigidDetermineCandidateNodesFromAdjacentComponents(dec, newRow, reducedId);
   if( SPQRnodeIsInvalid(nodes.first) && SPQRnodeIsInvalid(nodes.second))
   {
      newRow->remainsNetwork = FALSE;
      newRow->reducedMembers[reducedId].type = TYPE_NOT_NETWORK;
      return;
   }
   if( SPQRnodeIsInvalid(nodes.first) && SPQRnodeIsValid(nodes.second))
   {
      nodes.first = nodes.second;
      nodes.second = SPQR_INVALID_NODE;
   }

   spqr_node markerHead = findArcHead(dec, marker);
   spqr_node markerTail = findArcTail(dec, marker);
   if( findArcSign(dec, marker).reversed )
   {
      spqr_node temp = markerHead;
      markerHead = markerTail;
      markerTail = temp;
   }

   if( newRow->reducedMembers[reducedId].numCutArcs == 0 )
   {
      assert(SPQRnodeIsInvalid(nodes.second));//There must be at least two adjacent components
      if( nodes.first != markerHead && nodes.first != markerTail )
      {
         newRow->remainsNetwork = FALSE;
         newRow->reducedMembers[reducedId].type = TYPE_NOT_NETWORK;
         return;
      }
      SCIP_Bool reverse = previousOrientation.headSplit == ( nodes.first == markerTail );
      newRow->reducedMembers[reducedId].splitNode = nodes.first;
      newRow->reducedMembers[reducedId].otherIsSource = previousOrientation.otherIsSource;
      newRow->reducedMembers[reducedId].willBeReversed = reverse;
      newRow->reducedMembers[reducedId].type = TYPE_MERGED;
      return;
   }
   if( !SPQRnodeIsValid(newRow->reducedMembers[reducedId].splitNode))
   {
      assert(newRow->reducedMembers[reducedId].numCutArcs >
             0);//Checking for propagation only makes sense if there is at least one cut arc

      rigidFindStarNodes(dec, newRow, reducedId);
      if( !newRow->remainsNetwork )
      {
         return;
      }

      if( SPQRnodeIsInvalid(newRow->reducedMembers[reducedId].splitNode))
      {
         //not a star => attempt to find splittable nodes using articulation node algorithms
         //We save some work by telling the methods that only the marker nodes should be checked
         rigidGetSplittableArticulationPointsOnPath(dec, newRow, reducedId, nodes.first, nodes.second);
      }
      if( !newRow->remainsNetwork )
      {
         return;
      }
      if( SPQRnodeIsInvalid(newRow->reducedMembers[reducedId].splitNode))
      {
         newRow->remainsNetwork = FALSE;
         newRow->reducedMembers[reducedId].type = TYPE_NOT_NETWORK;
         return;
      }
   }
   spqr_node splitNode = determineAndColorSplitNode(dec, newRow, reducedId, nodes.first, nodes.second);

   if( SPQRnodeIsInvalid(splitNode))
   {
      newRow->remainsNetwork = FALSE;
      newRow->reducedMembers[reducedId].type = TYPE_NOT_NETWORK;
      return;
   }
   assert(splitNode == nodes.first || splitNode == nodes.second);
   assert(splitNode == markerHead || splitNode == markerTail);

   spqr_node otherNode = splitNode == markerHead ? markerTail : markerHead;
   SCIP_Bool headsMatch = previousOrientation.headSplit == ( splitNode == markerHead );

   COLOR_STATUS otherColor = newRow->nodeColors[otherNode];
   assert(otherColor == COLOR_SOURCE || otherColor == COLOR_SINK);
   SCIP_Bool otherIsSource = otherColor == COLOR_SOURCE;

   SCIP_Bool good = headsMatch == ( previousOrientation.otherIsSource == otherIsSource );

   if( !good )
   {
      newRow->remainsNetwork = FALSE;
      newRow->reducedMembers[reducedId].type = TYPE_NOT_NETWORK;
      return;
   }

   newRow->reducedMembers[reducedId].splitNode = splitNode;
   newRow->reducedMembers[reducedId].willBeReversed = !headsMatch;
   newRow->reducedMembers[reducedId].type = TYPE_MERGED;
}

static void determineSplitTypeNext(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id current,
   reduced_member_id next,
   SCIP_Bool currentIsParent
)
{
   spqr_member member = newRow->reducedMembers[next].member;
   SPQRMemberType type = getMemberType(dec, member);
   spqr_arc nextArc = currentIsParent ? markerToParent(dec, member) : markerOfParent(dec,
                                                                                     newRow->reducedMembers[current].member);
   spqr_arc currentArc = currentIsParent ? markerOfParent(dec, member) : markerToParent(dec,
                                                                                        newRow->reducedMembers[current].member);
   SplitOrientation orientation = getRelativeOrientation(dec, newRow, current, currentArc);
   assert(type == SPQR_MEMBERTYPE_PARALLEL || type == SPQR_MEMBERTYPE_RIGID || type == SPQR_MEMBERTYPE_SERIES);
   switch( type )
   {
      case SPQR_MEMBERTYPE_RIGID:
      {
         determineSplitTypeRigid(dec, newRow, next, member, nextArc, orientation);
         break;
      }
      case SPQR_MEMBERTYPE_PARALLEL:
      {
         determineSplitTypeParallel(dec, newRow, next, nextArc, orientation);
         break;
      }
      case SPQR_MEMBERTYPE_SERIES:
      {
         determineSplitTypeSeries(dec, newRow, next, nextArc, orientation);
         break;
      }

      case SPQR_MEMBERTYPE_LOOP:
      case SPQR_MEMBERTYPE_UNASSIGNED:
      default:
         newRow->remainsNetwork = FALSE;
         assert(FALSE);
   }
}

static void determineMergeableTypes(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id root
)
{
   assert(newRow->numReducedMembers <= newRow->memMergeTreeCallData);
   if( newRow->reducedMembers[root].numPropagatedChildren == newRow->reducedMembers[root].numChildren )
   {
      //Determine single component;
      if( newRow->reducedMembers[root].type == TYPE_UNDETERMINED )
      {
         spqr_member member = newRow->reducedMembers[root].member;
         switch( getMemberType(dec, member))
         {
            case SPQR_MEMBERTYPE_RIGID:
               determineSingleRowRigidType(dec, newRow, root);
               break;
            case SPQR_MEMBERTYPE_PARALLEL:
               determineSingleParallelType(dec, newRow, root);
               break;
            case SPQR_MEMBERTYPE_LOOP:
            case SPQR_MEMBERTYPE_SERIES:
               determineSingleSeriesType( newRow, root);
               break;
            case SPQR_MEMBERTYPE_UNASSIGNED:
            default:
               newRow->remainsNetwork = FALSE;
               assert(FALSE);
         }
      }
      return;
   }

   //go to a leaf. We need to start in a leaf to avoid the ambiguity of choosing an orientation
   //in members which have no cut arcs; otherwise, we might choose the wrong one
   reduced_member_id leaf = root;
   while( newRow->reducedMembers[leaf].numChildren != newRow->reducedMembers[leaf].numPropagatedChildren )
   {
      for( int i = 0; i < newRow->reducedMembers[leaf].numChildren; ++i )
      {
         children_idx idx = newRow->reducedMembers[leaf].firstChild + i;
         reduced_member_id child = newRow->childrenStorage[idx];
         if( newRow->reducedMembers[child].type != TYPE_PROPAGATED )
         {
            leaf = child;
            break;
         }
      }
   }
   determineSplitTypeFirstLeaf(dec, newRow, leaf);

   if( !newRow->remainsNetwork )
   {
      return;
   }
   reduced_member_id baseNode = leaf;
   reduced_member_id nextNode = newRow->reducedMembers[baseNode].parent;

   while( reducedMemberIsValid(nextNode))
   {
      //check this node
      determineSplitTypeNext(dec, newRow, baseNode, nextNode, FALSE);
      if( !newRow->remainsNetwork )
      {
         return;
      }
      //Add other nodes in the subtree
      //use a while loop to avoid recursion; we may get stack overflows for large graphs
      MergeTreeCallData * data = newRow->mergeTreeCallData;

      data[0].id = nextNode;
      data[0].currentChild = newRow->reducedMembers[nextNode].firstChild ;
      int depth = 0;
      while(depth >= 0){
         reduced_member_id id = data[depth].id;
         children_idx childidx = data[depth].currentChild;
         if(childidx == newRow->reducedMembers[id].firstChild + newRow->reducedMembers[id].numChildren){
            --depth;
            continue;
         }
         reduced_member_id currentchild = newRow->childrenStorage[childidx];
         data[depth].currentChild += 1;
         //skip this child if we already processed it or it is not merged
         if( currentchild == baseNode || newRow->reducedMembers[currentchild].type == TYPE_PROPAGATED){
            continue;
         }
         determineSplitTypeNext(dec,newRow,id,currentchild,TRUE);
         if(!newRow->remainsNetwork){
            return;
         }
         //recursively process the child
         depth += 1;
         assert(depth < newRow->memMergeTreeCallData);
         data[depth].id = currentchild;
         data[depth].currentChild = newRow->reducedMembers[currentchild].firstChild;
      }
      //Move up one layer
      baseNode = nextNode;
      nextNode = newRow->reducedMembers[nextNode].parent;
   }
}

static void cleanUpRowMemberInformation(SCIP_NETROWADD* newRow)
{
   //This loop is at the end as memberInformation is also used to assign the cut arcs during propagation
   //Clean up the memberInformation array
   for( int i = 0; i < newRow->numReducedMembers; ++i )
   {
      newRow->memberInformation[newRow->reducedMembers[i].member].reducedMember = INVALID_REDUCED_MEMBER;
   }
#ifndef NDEBUG
   for( int i = 0; i < newRow->memMemberInformation; ++i )
   {
      assert(reducedMemberIsInvalid(newRow->memberInformation[i].reducedMember));
   }
#endif
}

static SCIP_RETCODE rigidTransformArcIntoCycle(
   SCIP_NETMATDECDATA* dec,
   const spqr_member member,
   const spqr_arc arc,
   const SCIP_Bool reverseArcDirection,
   NewRowInformation* const newRowInformation
)
{
   //If a cycle already exists, just expand it with the new arc.
   spqr_member markerCycleMember = SPQR_INVALID_MEMBER;
   spqr_arc markerCycleArc = SPQR_INVALID_ARC;
   SCIP_Bool isParent = arc == markerToParent(dec, member);
   spqr_member adjacentMember = SPQR_INVALID_MEMBER;
   if( isParent )
   {
      adjacentMember = findMemberParent(dec, member);
      if( getMemberType(dec, adjacentMember) == SPQR_MEMBERTYPE_SERIES )
      {
         markerCycleMember = adjacentMember;
         markerCycleArc = markerOfParent(dec, member);
      }
   } else if( arcIsMarker(dec, arc))
   {
      adjacentMember = findArcChildMember(dec, arc);
      if( getMemberType(dec, adjacentMember) == SPQR_MEMBERTYPE_SERIES )
      {
         markerCycleMember = adjacentMember;
         markerCycleArc = markerToParent(dec, adjacentMember);
      }
   }
   if( markerCycleMember != SPQR_INVALID_MEMBER)
   {
      newRowInformation->member = markerCycleMember;
      if( arcIsReversedNonRigid(dec, markerCycleArc))
      {
         newRowInformation->reversed = reverseArcDirection;
      } else
      {
         newRowInformation->reversed = !reverseArcDirection;
      }

      return SCIP_OKAY;
   }

   //Otherwise, we create a new cycle
   spqr_member newCycle;
   SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &newCycle));
   //We would like to move the edge but unfortunately cannot do so without destroying the arc union-find datastructure.
   //Thus, we 'convert' the arc into a marker and add the new series

   spqr_arc duplicate = SPQR_INVALID_ARC;
   spqr_element element = arcGetElement(dec, arc);
   if( element != MARKER_COLUMN_ELEMENT && element != MARKER_ROW_ELEMENT)
   {
      if( SPQRelementIsColumn(element))
      {
         SCIP_CALL(createColumnArc(dec, newCycle, &duplicate, SPQRelementToColumn(element), TRUE));
      } else
      {
         SCIP_CALL(createRowArc(dec, newCycle, &duplicate, SPQRelementToRow(element), TRUE));
      }
   } else if( isParent )
   {
      //create parent marker
      SCIP_CALL(createParentMarker(dec, newCycle, arcIsTree(dec, arc), adjacentMember,
                                   markerOfParent(dec, member), &duplicate, TRUE));
   } else
   {
      //create child marker
      SCIP_CALL(createChildMarker(dec, newCycle, adjacentMember, arcIsTree(dec, arc), &duplicate, TRUE));
      dec->members[adjacentMember].parentMember = newCycle;
      dec->members[adjacentMember].markerOfParent = duplicate;
   }
   //Create the other marker edge
   spqr_arc cycleMarker = SPQR_INVALID_ARC;
   if( isParent )
   {
      SCIP_CALL(createChildMarker(dec, newCycle, member, !arcIsTree(dec, arc),
                                  &cycleMarker, FALSE));
   } else
   {
      SCIP_CALL(createParentMarker(dec, newCycle, !arcIsTree(dec, arc),
                                   member, arc, &cycleMarker, FALSE));
   }
   //Change the existing edge to a marker
   if( isParent )
   {
      assert(markerToParent(dec, member) == arc);
      dec->arcs[markerOfParent(dec, member)].childMember = newCycle;
      dec->members[member].parentMember = newCycle;
      dec->members[member].markerToParent = arc;
      dec->members[member].markerOfParent = cycleMarker;
      dec->arcs[arc].element = arcIsTree(dec, arc) ? MARKER_ROW_ELEMENT : MARKER_COLUMN_ELEMENT;
      dec->arcs[arc].childMember = SPQR_INVALID_MEMBER;

   } else
   {
      dec->arcs[arc].element = arcIsTree(dec, arc) ? MARKER_ROW_ELEMENT : MARKER_COLUMN_ELEMENT;
      dec->arcs[arc].childMember = newCycle;
   }
   newRowInformation->member = newCycle;
   newRowInformation->reversed = !reverseArcDirection;

   return SCIP_OKAY;
}

static SCIP_RETCODE transformSingleRigid(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const reduced_member_id reducedMember,
   const spqr_member member,
   NewRowInformation* const newRowInformation
)
{
   if( SPQRarcIsValid(newRow->reducedMembers[reducedMember].articulationArc))
   {
      spqr_arc arc = newRow->reducedMembers[reducedMember].articulationArc;
      //Cut arc is propagated into a cycle with new arc
      assert(newRow->reducedMembers[reducedMember].splitNode == findEffectiveArcHeadNoCompression(dec, arc) ||
             newRow->reducedMembers[reducedMember].splitNode == findEffectiveArcTailNoCompression(dec, arc));


      SCIP_Bool reversed;

      if( newRow->isArcCut[arc] )
      {
         reversed = ( newRow->reducedMembers[reducedMember].splitNode == findEffectiveArcHead(dec, arc)) ==
                    newRow->reducedMembers[reducedMember].otherIsSource;
      } else
      {
         reversed = ( newRow->reducedMembers[reducedMember].splitNode == findEffectiveArcHead(dec, arc)) !=
                    newRow->reducedMembers[reducedMember].otherIsSource;
      }

      SCIP_CALL(rigidTransformArcIntoCycle(dec, member, newRow->reducedMembers[reducedMember].articulationArc,
                                           reversed, newRowInformation));

      return SCIP_OKAY;
   }
   //Single splittable node
   assert(SPQRnodeIsValid(newRow->reducedMembers[reducedMember].splitNode));

   spqr_node splitNode = newRow->reducedMembers[reducedMember].splitNode;
   if( newRow->reducedMembers[reducedMember].allHaveCommonNode )
   {
      //Create a new node; move all cut arcs end of split node to it and add new arc between new node and split node
      spqr_node newNode = SPQR_INVALID_NODE;
      SCIP_CALL(createNode(dec, &newNode));

      cut_arc_id cutArcIdx = newRow->reducedMembers[reducedMember].firstCutArc;
      do
      {
         spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
         spqr_node arcHead = findArcHead(dec, cutArc);
         if( arcHead == splitNode )
         {
            changeArcHead(dec, cutArc, arcHead, newNode);
         } else
         {
            changeArcTail(dec, cutArc, findArcTail(dec, cutArc), newNode);
         }

         cutArcIdx = newRow->cutArcs[cutArcIdx].nextMember;
      } while( cutArcIsValid(cutArcIdx));

      newRowInformation->member = member;
      if( newRow->reducedMembers[reducedMember].otherIsSource )
      {
         newRowInformation->head = newNode;
         newRowInformation->tail = splitNode;
      } else
      {
         newRowInformation->head = splitNode;
         newRowInformation->tail = newNode;
      }
      newRowInformation->representative = findArcSign(dec,
                                                      newRow->cutArcs[newRow->reducedMembers[reducedMember].firstCutArc].arc).representative;
      newRowInformation->reversed = FALSE;

      return SCIP_OKAY;
   }
   //Articulation point was split (based on coloring)

   spqr_node newNode = SPQR_INVALID_NODE;
   SCIP_CALL(createNode(dec, &newNode));

   spqr_arc firstNodeArc = getFirstNodeArc(dec, splitNode);
   spqr_arc iterArc = firstNodeArc;

   do
   {
      SCIP_Bool isCut = newRow->isArcCut[iterArc];
      spqr_node otherHead = findArcHead(dec, iterArc);
      spqr_node otherTail = findArcTail(dec, iterArc);
      spqr_node otherEnd = otherHead == splitNode ? otherTail : otherHead;
      SCIP_Bool isMoveColor = newRow->nodeColors[otherEnd] == COLOR_SOURCE;
      spqr_arc nextArc = getNextNodeArc(dec, iterArc, splitNode);//Need to do this before we modify the arc :)

      SCIP_Bool changeArcEnd = isCut == isMoveColor;
      if( changeArcEnd )
      {
         if( otherHead == splitNode )
         {
            changeArcHead(dec, iterArc, otherHead, newNode);
         } else
         {
            changeArcTail(dec, iterArc, otherTail, newNode);
         }
      }
      newRow->nodeColors[otherEnd] = UNCOLORED;//Clean up

      //Ugly hack to make sure we can iterate neighbourhood whilst changing arc ends.
      spqr_arc previousArc = iterArc;
      iterArc = nextArc;
      if( iterArc == firstNodeArc )
      {
         break;
      }
      if( changeArcEnd && previousArc == firstNodeArc )
      {
         firstNodeArc = iterArc;
      }
   } while( TRUE ); /*lint !e506*/
   newRow->reducedMembers[reducedMember].coloredNode = SPQR_INVALID_NODE;

   newRowInformation->member = member;
   newRowInformation->head = newNode;
   newRowInformation->tail = splitNode;
   newRowInformation->representative = findArcSign(dec,
                                                   newRow->cutArcs[newRow->reducedMembers[reducedMember].firstCutArc].arc).representative;
   newRowInformation->reversed = FALSE;

   return SCIP_OKAY;
}


static SCIP_RETCODE splitParallelRowAddition(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const reduced_member_id reducedMember,
   const spqr_member member,
   NewRowInformation* newRowInfo
)
{
   assert(newRow->reducedMembers[reducedMember].numCutArcs > 0);

   int numCutArcs = newRow->reducedMembers[reducedMember].numCutArcs;
   int numParallelArcs = getNumMemberArcs(dec, member);

   SCIP_Bool createCutParallel = numCutArcs > 1;
   SCIP_Bool convertOriginalParallel = ( numCutArcs + 1 ) == numParallelArcs;

   //Do linear search to find non-marked arc
   spqr_arc treeArc = getFirstMemberArc(dec, member);
   do
   {
      if( arcIsTree(dec, treeArc))
      {
         break;
      }
      treeArc = getNextMemberArc(dec, treeArc);
   } while( treeArc != getFirstMemberArc(dec, member));
   assert(arcIsTree(dec, treeArc));

   SCIP_Bool treeReversed = arcIsReversedNonRigid(dec, treeArc);

   assert(!( !createCutParallel &&
             convertOriginalParallel ));//This can only happen if the parallel member is actually a loop, which means it is mislabeled
   if( createCutParallel )
   {
      if( convertOriginalParallel )
      {
         spqr_member adjacentMember = SPQR_INVALID_MEMBER;
         spqr_arc adjacentArc = SPQR_INVALID_ARC;
         if( treeArc == markerToParent(dec, member))
         {
            adjacentMember = findMemberParent(dec, member);
            adjacentArc = markerOfParent(dec, member);
         } else if( arcIsMarker(dec, treeArc))
         {
            adjacentMember = findArcChildMember(dec, treeArc);
            adjacentArc = markerToParent(dec, adjacentMember);
            assert(markerOfParent(dec, adjacentMember) == treeArc);
         }
         cut_arc_id firstCut = newRow->reducedMembers[reducedMember].firstCutArc;
         SCIP_Bool firstReversed =
            newRow->cutArcs[firstCut].arcReversed != arcIsReversedNonRigid(dec, newRow->cutArcs[firstCut].arc);

         if( SPQRmemberIsValid(adjacentMember) && getMemberType(dec, adjacentMember) == SPQR_MEMBERTYPE_SERIES )
         {
            newRowInfo->member = adjacentMember;
            if( arcIsReversedNonRigid(dec, treeArc) == arcIsReversedNonRigid(dec, adjacentArc))
            {
               newRowInfo->reversed = !firstReversed;
            } else
            {
               newRowInfo->reversed = firstReversed;
            }
            return SCIP_OKAY;
         }
         spqr_member cutMember = SPQR_INVALID_MEMBER;
         SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &cutMember));

         cut_arc_id cutArcIdx = newRow->reducedMembers[reducedMember].firstCutArc;
         assert(cutArcIsValid(cutArcIdx));
         SCIP_Bool parentCut = FALSE;

         while( cutArcIsValid(cutArcIdx))
         {
            spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
            cutArcIdx = newRow->cutArcs[cutArcIdx].nextMember;
            moveArcToNewMember(dec, cutArc, member, cutMember);
            if( cutArc == markerToParent(dec, member))
            {
               parentCut = TRUE;
            }
         }
         if( parentCut )
         {
            SCIP_CALL(createMarkerPair(dec, cutMember, member, TRUE, FALSE, TRUE));
         } else
         {
            SCIP_CALL(createMarkerPair(dec, member, cutMember, FALSE, TRUE, FALSE));
         }
         changeLoopToSeries(dec, member);
         newRowInfo->member = member;
         if( treeReversed )
         {
            newRowInfo->reversed = firstReversed == arcIsReversedNonRigid(dec, treeArc);
         } else
         {
            newRowInfo->reversed = firstReversed != arcIsReversedNonRigid(dec, treeArc);
         }
      } else
      {
         spqr_member cutMember = SPQR_INVALID_MEMBER;
         SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &cutMember));

         cut_arc_id cutArcIdx = newRow->reducedMembers[reducedMember].firstCutArc;
         assert(cutArcIsValid(cutArcIdx));
         SCIP_Bool parentCut = FALSE;

         while( cutArcIsValid(cutArcIdx))
         {
            spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
            cutArcIdx = newRow->cutArcs[cutArcIdx].nextMember;
            moveArcToNewMember(dec, cutArc, member, cutMember);
            if( cutArc == markerToParent(dec, member))
            {
               parentCut = TRUE;
            }
         }
         spqr_member newSeries;
         SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &newSeries));
         if( parentCut )
         {
            SCIP_CALL(createMarkerPair(dec, newSeries, member, TRUE, FALSE, FALSE));
            SCIP_CALL(createMarkerPair(dec, cutMember, newSeries, TRUE, FALSE, TRUE));
         } else
         {
            SCIP_CALL(createMarkerPair(dec, member, newSeries, FALSE, FALSE, FALSE));
            SCIP_CALL(createMarkerPair(dec, newSeries, cutMember, FALSE, TRUE, FALSE));
         }
         newRowInfo->member = newSeries;
         cut_arc_id firstCut = newRow->reducedMembers[reducedMember].firstCutArc;
         SCIP_Bool firstReversed =
            newRow->cutArcs[firstCut].arcReversed != arcIsReversedNonRigid(dec, newRow->cutArcs[firstCut].arc);

         newRowInfo->reversed = firstReversed;
      }

      return SCIP_OKAY;
   }

   assert(!createCutParallel && !convertOriginalParallel);
   assert(numCutArcs == 1);
#ifndef NDEBUG
   spqr_arc arc = newRow->cutArcs[newRow->reducedMembers[reducedMember].firstCutArc].arc;
   spqr_member adjacentMember = SPQR_INVALID_MEMBER;
   if( arc == markerToParent(dec, member))
   {
      adjacentMember = findMemberParent(dec, member);
   } else if( arcIsMarker(dec, arc))
   {
      adjacentMember = findArcChildMember(dec, arc);
   }
   if( SPQRmemberIsValid(adjacentMember))
   {
      assert(getMemberType(dec, adjacentMember) != SPQR_MEMBERTYPE_SERIES);
   }
#endif
   spqr_member newSeries;
   SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &newSeries));
   cut_arc_id cutArcIdx = newRow->reducedMembers[reducedMember].firstCutArc;
   assert(cutArcIsValid(cutArcIdx));
   SCIP_Bool parentCut = FALSE;

   while( cutArcIsValid(cutArcIdx))
   {
      spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
      cutArcIdx = newRow->cutArcs[cutArcIdx].nextMember;
      moveArcToNewMember(dec, cutArc, member, newSeries);
      if( cutArc == markerToParent(dec, member))
      {
         parentCut = TRUE;
      }
   }
   if( parentCut )
   {
      SCIP_CALL(createMarkerPair(dec, newSeries, member, TRUE, TRUE, FALSE));
   } else
   {
      SCIP_CALL(createMarkerPair(dec, member, newSeries, FALSE, FALSE, TRUE));
   }
   newRowInfo->member = newSeries;
   cut_arc_id cutArcId = newRow->reducedMembers[reducedMember].firstCutArc;
   newRowInfo->reversed =
      newRow->cutArcs[cutArcId].arcReversed == arcIsReversedNonRigid(dec, newRow->cutArcs[cutArcId].arc);

   return SCIP_OKAY;
}

static SCIP_RETCODE transformSingleParallel(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const reduced_member_id reducedMember,
   const spqr_member member,
   NewRowInformation* info
)
{
   SCIP_CALL(splitParallelRowAddition(dec, newRow, reducedMember, member, info));
   return SCIP_OKAY;
}

static SCIP_RETCODE splitSeriesMergingRowAddition(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   const reduced_member_id reducedMember,
   const spqr_member member,
   spqr_member* const mergingMember,
   SCIP_Bool* const isCut,
   spqr_arc* representativeEdge
)
{
   assert(getNumMemberArcs(dec, member) >= 3);
   *isCut = newRow->reducedMembers[reducedMember].numCutArcs > 0;

   if( getNumMemberArcs(dec, member) == 3 )
   {
      spqr_arc splitArc = newRow->reducedMembers[reducedMember].splitArc;
      spqr_arc otherArc = SPQR_INVALID_ARC;
      for( children_idx i = newRow->reducedMembers[reducedMember].firstChild;
           i <
           newRow->reducedMembers[reducedMember].firstChild + newRow->reducedMembers[reducedMember].numChildren; ++i )
      {
         reduced_member_id child = newRow->childrenStorage[i];
         if( newRow->reducedMembers[child].type == TYPE_MERGED )
         {
            spqr_arc testArc = markerOfParent(dec, findMember(dec, newRow->reducedMembers[child].member));
            if( testArc != splitArc )
            {
               otherArc = testArc;
               break;
            }
         }
      }
      if( SPQRarcIsInvalid(otherArc))
      {
#ifndef NDEBUG
         reduced_member_id parent = newRow->reducedMembers[reducedMember].parent;
#endif
         assert(newRow->reducedMembers[parent].type == TYPE_MERGED ||
                newRow->reducedMembers[parent].type == TYPE_PROPAGATED);
         assert(reducedMemberIsValid(parent) && newRow->reducedMembers[parent].type == TYPE_MERGED);
         otherArc = markerToParent(dec, member);
      }
      spqr_arc firstArc = getFirstMemberArc(dec, member);
      spqr_arc arc = firstArc;
      do
      {
         if( arc != splitArc && arc != otherArc )
         {
            *representativeEdge = arc;
            break;
         }
         arc = getNextMemberArc(dec, arc);
      } while( arc != firstArc );
      *mergingMember = member;
      return SCIP_OKAY;
   }
   //Split off the relevant part of the series member
   spqr_member mergingSeries = SPQR_INVALID_MEMBER;
   SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_SERIES, &mergingSeries));
   //Move all marker arcs which point to another component in the reduced decomposition to the new member
   //This should be exactly 2, as with 3 the result is not network anymore
   //move all mergeable children and parent arcs to the mergingMember

   SCIP_Bool coTreeToMergingMember = FALSE;
   SCIP_Bool parentToMergingMember = FALSE;
   for( children_idx i = newRow->reducedMembers[reducedMember].firstChild;
        i < newRow->reducedMembers[reducedMember].firstChild + newRow->reducedMembers[reducedMember].numChildren; ++i )
   {
      reduced_member_id child = newRow->childrenStorage[i];
      assert(
         newRow->reducedMembers[child].type == TYPE_MERGED || newRow->reducedMembers[child].type == TYPE_PROPAGATED);
      if( newRow->reducedMembers[child].type == TYPE_MERGED )
      {
         spqr_arc moveArc = markerOfParent(dec, findMember(dec, newRow->reducedMembers[child].member));
         moveArcToNewMember(dec, moveArc, member, mergingSeries);
         if( !arcIsTree(dec, moveArc))
         {
            coTreeToMergingMember = TRUE;
         }
      }
   }

   reduced_member_id parent = newRow->reducedMembers[reducedMember].parent;
   assert(reducedMemberIsInvalid(parent) || ( newRow->reducedMembers[parent].type == TYPE_MERGED ||
                                              newRow->reducedMembers[parent].type == TYPE_PROPAGATED ));

   if( reducedMemberIsValid(parent) &&
       newRow->reducedMembers[parent].type == TYPE_MERGED )
   {
      spqr_arc moveArc = markerToParent(dec, member);
      moveArcToNewMember(dec, moveArc, member, mergingSeries);
      parentToMergingMember = TRUE;
      if( !arcIsTree(dec, moveArc))
      {
         coTreeToMergingMember = TRUE;
      }
   }
   spqr_arc ignoreArc = SPQR_INVALID_ARC;
   if( parentToMergingMember )
   {
      SCIP_CALL(createMarkerPairWithReferences(dec, mergingSeries, member, coTreeToMergingMember, TRUE, FALSE,
                                               representativeEdge, &ignoreArc));
   } else
   {
      SCIP_CALL(
         createMarkerPairWithReferences(dec, member, mergingSeries, !coTreeToMergingMember, FALSE, TRUE, &ignoreArc,
                                        representativeEdge));
   }

   *mergingMember = mergingSeries;
   assert(getNumMemberArcs(dec, mergingSeries) == 3);
   assert(getNumMemberArcs(dec, member) >= 3);
   return SCIP_OKAY;
}

static SCIP_RETCODE splitParallelMerging(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id reducedMember,
   spqr_member member,
   spqr_member* const pMergeMember,
   spqr_arc* const cutRepresentative
)
{
   //When merging, we cannot have propagated members;
   assert(newRow->reducedMembers[reducedMember].numCutArcs < ( getNumMemberArcs(dec, member) - 1 ));

   int numMergeableAdjacent =
      newRow->reducedMembers[reducedMember].numChildren - newRow->reducedMembers[reducedMember].numPropagatedChildren;
   if( reducedMemberIsValid(newRow->reducedMembers[reducedMember].parent) &&
       newRow->reducedMembers[newRow->reducedMembers[reducedMember].parent].type == TYPE_MERGED )
   {
      numMergeableAdjacent++;
   }

   int numCutArcs = newRow->reducedMembers[reducedMember].numCutArcs;
   //All arcs which are not in the mergeable decomposition or cut
   int numBaseSplitAwayArcs = getNumMemberArcs(dec, member) - numMergeableAdjacent - numCutArcs;

   SCIP_Bool createCutParallel = numCutArcs > 1;
   SCIP_Bool keepOriginalParallel = numBaseSplitAwayArcs <= 1;

   spqr_member cutMember = SPQR_INVALID_MEMBER;
   //The below cases can probably be aggregated in some way, but for now we first focus on the correct logic
   if( createCutParallel && keepOriginalParallel )
   {
      SCIP_Bool parentCut = FALSE;
      SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &cutMember));

      cut_arc_id cutArcIdx = newRow->reducedMembers[reducedMember].firstCutArc;
      assert(cutArcIsValid(cutArcIdx));

      while( cutArcIsValid(cutArcIdx))
      {
         spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
         cutArcIdx = newRow->cutArcs[cutArcIdx].nextMember;
         moveArcToNewMember(dec, cutArc, member, cutMember);
         if( cutArc == markerToParent(dec, member))
         {
            parentCut = TRUE;
         }
      }
      spqr_arc ignoreArc = SPQR_INVALID_ARC;
      if( parentCut )
      {
         SCIP_CALL(
            createMarkerPairWithReferences(dec, cutMember, member, TRUE, FALSE, FALSE, &ignoreArc, cutRepresentative));
      } else
      {
         SCIP_CALL(
            createMarkerPairWithReferences(dec, member, cutMember, FALSE, FALSE, FALSE, cutRepresentative, &ignoreArc));
      }

      *pMergeMember = member;
   } else if( createCutParallel )
   {
      assert(!keepOriginalParallel);

      SCIP_Bool parentCut = FALSE;
      SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &cutMember));

      cut_arc_id cutArcIdx = newRow->reducedMembers[reducedMember].firstCutArc;
      assert(cutArcIsValid(cutArcIdx));

      while( cutArcIsValid(cutArcIdx))
      {
         spqr_arc cutArc = newRow->cutArcs[cutArcIdx].arc;
         cutArcIdx = newRow->cutArcs[cutArcIdx].nextMember;
         moveArcToNewMember(dec, cutArc, member, cutMember);
         if( cutArc == markerToParent(dec, member))
         {
            parentCut = TRUE;
         }
      }
      spqr_arc ignoreArc = SPQR_INVALID_ARC;
      if( parentCut )
      {
         SCIP_CALL(
            createMarkerPairWithReferences(dec, cutMember, member, TRUE, FALSE, FALSE, &ignoreArc, cutRepresentative));
      } else
      {
         SCIP_CALL(
            createMarkerPairWithReferences(dec, member, cutMember, FALSE, FALSE, FALSE, cutRepresentative, &ignoreArc));
      }


      spqr_arc noCutRepresentative = SPQR_INVALID_ARC;
      spqr_member mergingMember = member;
      SCIP_Bool parentToMergingMember = FALSE;
      SCIP_Bool treeToMergingMember = FALSE;
      SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &mergingMember));
      //move all mergeable children and parent arcs to the mergingMember
      for( children_idx i = newRow->reducedMembers[reducedMember].firstChild;
           i < newRow->reducedMembers[reducedMember].firstChild +
               newRow->reducedMembers[reducedMember].numChildren;
           ++i )
      {
         reduced_member_id child = newRow->childrenStorage[i];
         assert(
            newRow->reducedMembers[child].type == TYPE_MERGED || newRow->reducedMembers[child].type == TYPE_PROPAGATED);
         if( newRow->reducedMembers[child].type == TYPE_MERGED )
         {
            spqr_arc moveArc = markerOfParent(dec, findMember(dec, newRow->reducedMembers[child].member));
            moveArcToNewMember(dec, moveArc, member, mergingMember);
            if( arcIsTree(dec, moveArc))
            {
               treeToMergingMember = TRUE;
            }
         }
      }
      reduced_member_id parent = newRow->reducedMembers[reducedMember].parent;
      if( reducedMemberIsValid(parent) &&
          newRow->reducedMembers[parent].type == TYPE_MERGED )
      {
         spqr_arc moveArc = markerToParent(dec, member);
         moveArcToNewMember(dec, moveArc, member, mergingMember);
         parentToMergingMember = TRUE;
         if( arcIsTree(dec, moveArc))
         {
            treeToMergingMember = TRUE;
         }
      }
      //If there is only one cut arc, we also move it.
      if( SPQRarcIsValid(*cutRepresentative))
      {
         if( *cutRepresentative == markerToParent(dec, member))
         {
            parentToMergingMember = TRUE;
         }
         moveArcToNewMember(dec, *cutRepresentative, member, mergingMember);
      }
      spqr_arc ignoreArgument = SPQR_INVALID_ARC;
      if( parentToMergingMember || parentCut )
      {
         SCIP_CALL(createMarkerPairWithReferences(dec, mergingMember, member, !treeToMergingMember, FALSE, FALSE,
                                                  &ignoreArgument, &noCutRepresentative));
      } else
      {
         SCIP_CALL(createMarkerPairWithReferences(dec, member, mergingMember, treeToMergingMember, FALSE, FALSE,
                                                  &noCutRepresentative, &ignoreArgument));
      }
      *pMergeMember = mergingMember;
   } else if( keepOriginalParallel )
   {
      assert(!createCutParallel);
      if( cutArcIsValid(newRow->reducedMembers[reducedMember].firstCutArc))
      {
         *cutRepresentative = newRow->cutArcs[newRow->reducedMembers[reducedMember].firstCutArc].arc;
      }
      *pMergeMember = member;

   } else
   {
      assert(!keepOriginalParallel && !createCutParallel);

      if( cutArcIsValid(newRow->reducedMembers[reducedMember].firstCutArc))
      {
         *cutRepresentative = newRow->cutArcs[newRow->reducedMembers[reducedMember].firstCutArc].arc;
      }

      spqr_arc noCutRepresentative = SPQR_INVALID_ARC;
      spqr_member mergingMember = member;
      SCIP_Bool parentToMergingMember = FALSE;
      SCIP_Bool treeToMergingMember = FALSE;
      SCIP_CALL(createMember(dec, SPQR_MEMBERTYPE_PARALLEL, &mergingMember));
      //move all mergeable children and parent arcs to the mergingMember
      for( children_idx i = newRow->reducedMembers[reducedMember].firstChild;
           i < newRow->reducedMembers[reducedMember].firstChild +
               newRow->reducedMembers[reducedMember].numChildren;
           ++i )
      {
         reduced_member_id child = newRow->childrenStorage[i];
         assert(
            newRow->reducedMembers[child].type == TYPE_MERGED || newRow->reducedMembers[child].type == TYPE_PROPAGATED);
         if( newRow->reducedMembers[child].type == TYPE_MERGED )
         {
            spqr_arc moveArc = markerOfParent(dec, findMember(dec, newRow->reducedMembers[child].member));
            moveArcToNewMember(dec, moveArc, member, mergingMember);
            if( arcIsTree(dec, moveArc))
            {
               treeToMergingMember = TRUE;
            }
         }
      }
      reduced_member_id parent = newRow->reducedMembers[reducedMember].parent;
      if( reducedMemberIsValid(parent) &&
          newRow->reducedMembers[parent].type == TYPE_MERGED )
      {
         spqr_arc moveArc = markerToParent(dec, member);
         moveArcToNewMember(dec, moveArc, member, mergingMember);
         parentToMergingMember = TRUE;
         if( arcIsTree(dec, moveArc))
         {
            treeToMergingMember = TRUE;
         }
      }
      //If there is only one cut arc, we also move it.
      if( SPQRarcIsValid(*cutRepresentative))
      {
         if( *cutRepresentative == markerToParent(dec, member))
         {
            parentToMergingMember = TRUE;
         }
         moveArcToNewMember(dec, *cutRepresentative, member, mergingMember);
      }
      spqr_arc ignoreArgument = SPQR_INVALID_ARC;
      if( parentToMergingMember )
      {
         SCIP_CALL(createMarkerPairWithReferences(dec, mergingMember, member, !treeToMergingMember, FALSE, FALSE,
                                                  &ignoreArgument, &noCutRepresentative));
      } else
      {
         SCIP_CALL(createMarkerPairWithReferences(dec, member, mergingMember, treeToMergingMember, FALSE, FALSE,
                                                  &noCutRepresentative, &ignoreArgument));
      }
      *pMergeMember = mergingMember;
   }
   return SCIP_OKAY;
}

static SCIP_RETCODE splitFirstLeaf(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id leaf,
   NewRowInformation* const newRowInformation
)
{
   assert(cutArcIsValid(newRow->reducedMembers[leaf].firstCutArc));
   assert(newRow->reducedMembers[leaf].numCutArcs > 0);
   spqr_member member = newRow->reducedMembers[leaf].member;
   if( getMemberType(dec, member) == SPQR_MEMBERTYPE_PARALLEL )
   {
      spqr_member mergeMember = SPQR_INVALID_MEMBER;
      spqr_arc cutRepresentative = SPQR_INVALID_ARC;
      SCIP_CALL(splitParallelMerging(dec, newRow, leaf, member, &mergeMember, &cutRepresentative));
      newRow->reducedMembers[leaf].member = mergeMember;

      spqr_node firstNode = SPQR_INVALID_NODE;
      spqr_node secondNode = SPQR_INVALID_NODE;
      spqr_node thirdNode = SPQR_INVALID_NODE;
      SCIP_CALL(createNode(dec, &firstNode));
      SCIP_CALL(createNode(dec, &secondNode));
      SCIP_CALL(createNode(dec, &thirdNode));

      spqr_arc splitArc = newRow->reducedMembers[leaf].splitArc;
      assert(findArcMemberNoCompression(dec, splitArc) == mergeMember);
      SCIP_Bool splitArcReversed = arcIsReversedNonRigid(dec, splitArc);
      SCIP_Bool splitHead = newRow->reducedMembers[leaf].splitHead;
      spqr_node splitArcHead = splitArcReversed ? secondNode : firstNode;
      spqr_node splitArcTail = splitArcReversed ? firstNode : secondNode;
      spqr_node splitNode = splitHead ? splitArcHead : splitArcTail;
      spqr_node otherNode = splitHead ? splitArcTail : splitArcHead;

      spqr_arc first_arc = getFirstMemberArc(dec, mergeMember);
      spqr_arc arc = first_arc;

      do
      {
         if( arc != cutRepresentative )
         {
            if( arcIsReversedNonRigid(dec, arc))
            {
               setArcHeadAndTail(dec, arc, secondNode, firstNode);
            } else
            {
               setArcHeadAndTail(dec, arc, firstNode, secondNode);
            }
         } else
         {
            if(( arcIsReversedNonRigid(dec, arc) == splitArcReversed ) == splitHead )
            {
               setArcHeadAndTail(dec, arc, thirdNode, otherNode);
            } else
            {
               setArcHeadAndTail(dec, arc, otherNode, thirdNode);
            }
         }
         arcSetReversed(dec, arc, FALSE);
         if( arc == splitArc )
         {
            arcSetRepresentative(dec, arc, SPQR_INVALID_ARC);
         } else
         {
            arcSetRepresentative(dec, arc, splitArc);
         }
         arc = getNextMemberArc(dec, arc);
      } while( arc != first_arc );


      updateMemberType(dec, mergeMember, SPQR_MEMBERTYPE_RIGID);
      newRowInformation->member = mergeMember;
      if( newRow->reducedMembers[leaf].otherIsSource )
      {
         newRowInformation->head = thirdNode;
         newRowInformation->tail = splitNode;
      } else
      {
         newRowInformation->head = splitNode;
         newRowInformation->tail = thirdNode;
      }

      newRowInformation->reversed = FALSE;
      newRowInformation->representative = splitArc;

      return SCIP_OKAY;
   }
   assert(getMemberType(dec, member) == SPQR_MEMBERTYPE_RIGID);

   spqr_node newNode = SPQR_INVALID_NODE;//Sink node
   SCIP_CALL(createNode(dec, &newNode));

   spqr_node splitNode = newRow->reducedMembers[leaf].splitNode;

   spqr_arc firstNodeArc = getFirstNodeArc(dec, splitNode);
   spqr_arc iterArc = firstNodeArc;
   do
   {
      spqr_node otherHead = findArcHead(dec, iterArc);
      spqr_node otherTail = findArcTail(dec, iterArc);
      spqr_node otherEnd = otherHead == splitNode ? otherTail : otherHead;
      spqr_arc nextArc = getNextNodeArc(dec, iterArc, splitNode);//Need to do this before we modify the arc

      SCIP_Bool isCut = newRow->isArcCut[iterArc];
      SCIP_Bool isMoveColor = newRow->nodeColors[otherEnd] == COLOR_SOURCE;
      SCIP_Bool changeArcEnd = isCut == isMoveColor;
      if( changeArcEnd )
      {
         if( otherHead == splitNode )
         {
            changeArcHead(dec, iterArc, otherHead, newNode);
         } else
         {
            changeArcTail(dec, iterArc, otherTail, newNode);
         }
      }
      //Ugly hack to make sure we can iterate neighbourhood whilst changing arc ends.
      newRow->nodeColors[otherEnd] = UNCOLORED;
      spqr_arc previousArc = iterArc;
      iterArc = nextArc;
      if( iterArc == firstNodeArc )
      {
         break;
      }
      if( changeArcEnd && previousArc == firstNodeArc )
      {
         firstNodeArc = iterArc;
      }
   } while( TRUE ); /*lint !e506*/
   newRowInformation->head = newNode;
   newRowInformation->tail = splitNode;
   newRowInformation->member = member;
   newRowInformation->reversed = FALSE;
   newRowInformation->representative = findArcSign(dec, iterArc).representative;

   return SCIP_OKAY;
}

static SCIP_RETCODE mergeSplitMemberIntoParent(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   spqr_member member,
   spqr_member parent,
   spqr_arc parentToChild,
   spqr_arc childToParent,
   SCIP_Bool headToHead,
   spqr_node parentNode,
   spqr_node childNode,
   spqr_member* mergedMember,
   spqr_node* arcNodeOne,
   spqr_node* arcNodeTwo,
   spqr_node* thirdNode
)
{
   assert(dec);
   assert(SPQRmemberIsValid(member));
   assert(memberIsRepresentative(dec, member));
   assert(SPQRmemberIsValid(parent));
   assert(memberIsRepresentative(dec, parent));
   assert(findMemberParentNoCompression(dec, member) == parent);
   assert(markerOfParent(dec, member) == parentToChild);
   assert(markerToParent(dec, member) == childToParent);

   removeArcFromMemberArcList(dec, parentToChild, parent);
   removeArcFromMemberArcList(dec, childToParent, member);

   int parentTail = findEffectiveArcTail(dec, parentToChild);
   int parentHead = findEffectiveArcHead(dec, parentToChild);
   int childTail = findEffectiveArcTail(dec, childToParent);
   int childHead = findEffectiveArcHead(dec, childToParent);
   spqr_node parentArcNodes[2] = { parentTail, parentHead };
   spqr_node childArcNodes[2] = { childTail, childHead };

   clearArcHeadAndTail(dec, parentToChild);
   clearArcHeadAndTail(dec, childToParent);

   spqr_node first = childArcNodes[headToHead ? 0 : 1];
   spqr_node second = childArcNodes[headToHead ? 1 : 0];
   {
      spqr_node newNode = mergeNodes(dec, parentArcNodes[0], first);
      spqr_node toRemoveFrom = newNode == first ? parentArcNodes[0] : first;
      mergeNodeArcList(dec, newNode, toRemoveFrom);
      *arcNodeOne = newNode;
      newRow->nodeColors[toRemoveFrom] = UNCOLORED;
   }
   {
      spqr_node newNode = mergeNodes(dec, parentArcNodes[1], second);
      spqr_node toRemoveFrom = newNode == second ? parentArcNodes[1] : second;
      mergeNodeArcList(dec, newNode, toRemoveFrom);
      *arcNodeTwo = newNode;
      newRow->nodeColors[toRemoveFrom] = UNCOLORED;
   }
   {
      spqr_node newNode = mergeNodes(dec, parentNode, childNode);
      spqr_node toRemoveFrom = newNode == childNode ? parentNode : childNode;
      mergeNodeArcList(dec, newNode, toRemoveFrom);
      *thirdNode = newNode;
      newRow->nodeColors[toRemoveFrom] = UNCOLORED;
   }

   spqr_member newMember = mergeMembers(dec, member, parent);
   spqr_member toRemoveFrom = newMember == member ? parent : member;
   mergeMemberArcList(dec, newMember, toRemoveFrom);
   if( toRemoveFrom == parent )
   {
      updateMemberParentInformation(dec, newMember, toRemoveFrom);
   }
   updateMemberType(dec, newMember, SPQR_MEMBERTYPE_RIGID);
   *mergedMember = newMember;
   return SCIP_OKAY;
}

static SCIP_RETCODE splitAndMergeSeries(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id smallMember,
   SCIP_Bool largeIsParent,
   NewRowInformation* const newRowInformation,
   spqr_member member
)
{
   SCIP_Bool isCut = FALSE;
   spqr_member mergingMember = SPQR_INVALID_MEMBER;
   spqr_arc nonVirtualArc = SPQR_INVALID_ARC;
   SCIP_CALL(splitSeriesMergingRowAddition(dec, newRow, smallMember, member, &mergingMember, &isCut, &nonVirtualArc));
   assert(getNumMemberArcs(dec, mergingMember) == 3);

   //create the split series. There's two possible configurations, based on whether it contains a cut edge or not
   spqr_node a = SPQR_INVALID_NODE;
   spqr_node b = SPQR_INVALID_NODE;
   spqr_node c = SPQR_INVALID_NODE;
   spqr_node d = SPQR_INVALID_NODE;
   SCIP_CALL(createNode(dec, &a));
   SCIP_CALL(createNode(dec, &b));
   SCIP_CALL(createNode(dec, &c));
   SCIP_CALL(createNode(dec, &d));

   spqr_arc splitArc = newRow->reducedMembers[smallMember].splitArc;

   {
      SCIP_Bool splitHead = newRow->reducedMembers[smallMember].splitHead;
      spqr_arc firstArc = getFirstMemberArc(dec, mergingMember);
      spqr_arc arc = firstArc;
      SCIP_Bool splitReversed = arcIsReversedNonRigid(dec, splitArc);
      do
      {
         if( arc == splitArc )
         {
            if( splitHead )
            {
               setArcHeadAndTail(dec, splitArc, b, a);
            } else
            {
               setArcHeadAndTail(dec, splitArc, a, b);
            }
         } else if( arc == nonVirtualArc )
         {
            if(( arcIsReversedNonRigid(dec, arc) == splitReversed ) == splitHead )
            {
               setArcHeadAndTail(dec, arc, a, d);
            } else
            {
               setArcHeadAndTail(dec, arc, d, a);
            }
         } else
         {
            spqr_node otherNode = cutArcIsValid(newRow->reducedMembers[smallMember].firstCutArc) ? c : b;
            if(( arcIsReversedNonRigid(dec, arc) == splitReversed ) == splitHead )
            {
               setArcHeadAndTail(dec, arc, d, otherNode);
            } else
            {
               setArcHeadAndTail(dec, arc, otherNode, d);
            }
         }
         arcSetReversed(dec, arc, FALSE);
         arcSetRepresentative(dec, arc, splitArc);
         arc = getNextMemberArc(dec, arc);
      } while( arc != firstArc );
      arcSetRepresentative(dec, splitArc, SPQR_INVALID_ARC);
   }

   spqr_member otherMember = newRowInformation->member;
   spqr_arc otherMarker = largeIsParent ? markerOfParent(dec, mergingMember) : markerToParent(dec, otherMember);

   assert(nodeIsRepresentative(dec, newRowInformation->tail));
   assert(nodeIsRepresentative(dec, newRowInformation->head));
   spqr_node splitNode = newRow->reducedMembers[smallMember].splitHead ? findEffectiveArcHead(dec, otherMarker)
                                                                       : findEffectiveArcTail(dec, otherMarker);

   spqr_node otherNode = splitNode == newRowInformation->head ? newRowInformation->tail : newRowInformation->head;
   assert(splitNode == newRowInformation->head || splitNode == newRowInformation->tail);
   newRowInformation->representative = mergeArcSigns(dec, newRowInformation->representative, splitArc, FALSE);

   spqr_member mergedMember = SPQR_INVALID_MEMBER;
   spqr_node arcNodeOne;
   spqr_node arcNodeTwo;
   spqr_node thirdNode;
   if( largeIsParent )
   {
      SCIP_CALL(mergeSplitMemberIntoParent(dec, newRow, mergingMember, otherMember, otherMarker, splitArc, TRUE,
                                           otherNode, c, &mergedMember,
                                           &arcNodeOne,
                                           &arcNodeTwo,
                                           &thirdNode));
   } else
   {
      SCIP_CALL(mergeSplitMemberIntoParent(dec, newRow, otherMember, mergingMember, splitArc, otherMarker, TRUE,
                                           c, otherNode, &mergedMember,
                                           &arcNodeOne,
                                           &arcNodeTwo,
                                           &thirdNode));
   }
   newRow->reducedMembers[smallMember].member = mergedMember;

   newRowInformation->member = mergedMember;
   SCIP_Bool splitIsReferenceHead = newRow->reducedMembers[smallMember].splitHead;
   SCIP_Bool splitIsNewRowHead = splitNode == newRowInformation->head;
   if( !splitIsReferenceHead && !splitIsNewRowHead )
   {
      newRowInformation->head = thirdNode;
      newRowInformation->tail = arcNodeOne;
   } else if( !splitIsReferenceHead && splitIsNewRowHead )
   {
      newRowInformation->head = arcNodeOne;
      newRowInformation->tail = thirdNode;
   } else if( splitIsReferenceHead && !splitIsNewRowHead )
   {
      newRowInformation->head = thirdNode;
      newRowInformation->tail = arcNodeTwo;
   } else if( splitIsReferenceHead && splitIsNewRowHead )
   {
      newRowInformation->head = arcNodeTwo;
      newRowInformation->tail = thirdNode;
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE splitAndMergeParallel(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id smallMember,
   SCIP_Bool largeIsParent,
   NewRowInformation* const newRowInformation,
   spqr_member member
)
{
   spqr_member mergeMember = SPQR_INVALID_MEMBER;
   spqr_arc cutRepresentative = SPQR_INVALID_ARC;
   SCIP_CALL(splitParallelMerging(dec, newRow, smallMember, member, &mergeMember, &cutRepresentative));
   newRow->reducedMembers[smallMember].member = mergeMember;

   spqr_node firstNode = SPQR_INVALID_NODE;
   spqr_node secondNode = SPQR_INVALID_NODE;
   spqr_node thirdNode = SPQR_INVALID_NODE;
   SCIP_CALL(createNode(dec, &firstNode));
   SCIP_CALL(createNode(dec, &secondNode));
   SCIP_CALL(createNode(dec, &thirdNode));

   spqr_arc splitArc = newRow->reducedMembers[smallMember].splitArc;
   assert(findArcMemberNoCompression(dec, splitArc) == mergeMember);
   SCIP_Bool splitArcReversed = arcIsReversedNonRigid(dec, splitArc);
   SCIP_Bool splitHead = newRow->reducedMembers[smallMember].splitHead;

   spqr_node splitArcHead = splitArcReversed ? secondNode : firstNode;
   spqr_node splitArcTail = splitArcReversed ? firstNode : secondNode;
   spqr_node otherNode = splitHead ? splitArcTail : splitArcHead;

   spqr_arc first_arc = getFirstMemberArc(dec, mergeMember);
   spqr_arc arc = first_arc;

   do
   {
      if( arc != cutRepresentative )
      {
         if( arcIsReversedNonRigid(dec, arc))
         {
            setArcHeadAndTail(dec, arc, secondNode, firstNode);
         } else
         {
            setArcHeadAndTail(dec, arc, firstNode, secondNode);
         }
      } else
      {
         if(( arcIsReversedNonRigid(dec, arc) == splitArcReversed ) == splitHead )
         {
            setArcHeadAndTail(dec, arc, thirdNode, otherNode);
         } else
         {
            setArcHeadAndTail(dec, arc, otherNode, thirdNode);
         }
      }
      arcSetReversed(dec, arc, FALSE);
      arcSetRepresentative(dec, arc, splitArc);
      arc = getNextMemberArc(dec, arc);
   } while( arc != first_arc );
   arcSetRepresentative(dec, splitArc, SPQR_INVALID_ARC);

   spqr_member otherMember = newRowInformation->member;
   spqr_arc otherMarker = largeIsParent ? markerOfParent(dec, mergeMember) : markerToParent(dec, otherMember);

   assert(nodeIsRepresentative(dec, newRowInformation->tail));
   assert(nodeIsRepresentative(dec, newRowInformation->head));
   spqr_node largeSplitNode = newRow->reducedMembers[smallMember].splitHead ? findEffectiveArcHead(dec, otherMarker)
                                                                            : findEffectiveArcTail(dec, otherMarker);

   spqr_node largeOtherNode =
      largeSplitNode == newRowInformation->head ? newRowInformation->tail : newRowInformation->head;
   assert(largeSplitNode == newRowInformation->head || largeSplitNode == newRowInformation->tail);

   newRowInformation->representative = mergeArcSigns(dec, newRowInformation->representative, splitArc, FALSE);

   spqr_member mergedMember = SPQR_INVALID_MEMBER;
   spqr_node arcNodeOne;
   spqr_node arcNodeTwo;
   spqr_node mergeNodeThree;
   if( largeIsParent )
   {
      SCIP_CALL(mergeSplitMemberIntoParent(dec, newRow, mergeMember, otherMember, otherMarker, splitArc, TRUE,
                                           largeOtherNode, thirdNode, &mergedMember,
                                           &arcNodeOne,
                                           &arcNodeTwo,
                                           &mergeNodeThree));
   } else
   {
      SCIP_CALL(mergeSplitMemberIntoParent(dec, newRow, otherMember, mergeMember, splitArc, otherMarker, TRUE,
                                           thirdNode, largeOtherNode, &mergedMember,
                                           &arcNodeOne,
                                           &arcNodeTwo,
                                           &mergeNodeThree));
   }


   newRowInformation->member = mergedMember;

   SCIP_Bool splitIsReferenceHead = newRow->reducedMembers[smallMember].splitHead;
   SCIP_Bool splitIsNewRowHead = largeSplitNode == newRowInformation->head;
   if( !splitIsReferenceHead && !splitIsNewRowHead )
   {
      newRowInformation->head = mergeNodeThree;
      newRowInformation->tail = arcNodeOne;
   } else if( !splitIsReferenceHead && splitIsNewRowHead )
   {
      newRowInformation->head = arcNodeOne;
      newRowInformation->tail = mergeNodeThree;
   } else if( splitIsReferenceHead && !splitIsNewRowHead )
   {
      newRowInformation->head = mergeNodeThree;
      newRowInformation->tail = arcNodeTwo;
   } else if( splitIsReferenceHead && splitIsNewRowHead )
   {
      newRowInformation->head = arcNodeTwo;
      newRowInformation->tail = mergeNodeThree;
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE splitAndMergeRigid(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id smallMember,
   SCIP_Bool largeIsParent,
   NewRowInformation* const newRowInformation,
   spqr_member member
)
{

   spqr_node newNode = SPQR_INVALID_NODE;//Sink node
   SCIP_CALL(createNode(dec, &newNode));

   spqr_member smallMemberMember = member;
   spqr_member largeMemberMember = newRowInformation->member;

   spqr_arc smallMarker = largeIsParent ? markerToParent(dec, smallMemberMember) : markerOfParent(dec,
                                                                                                  largeMemberMember);
   spqr_arc largeMarker = largeIsParent ? markerOfParent(dec, smallMemberMember) : markerToParent(dec,
                                                                                                  largeMemberMember);

   spqr_node splitNode = newRow->reducedMembers[smallMember].splitNode;
   spqr_node smallOtherNode = newNode;

   if( newRow->reducedMembers[smallMember].numCutArcs != 0 )
   {
      spqr_arc firstNodeArc = getFirstNodeArc(dec, splitNode);
      spqr_arc iterArc = firstNodeArc;
      do
      {
         spqr_node otherHead = findArcHead(dec, iterArc);
         spqr_node otherTail = findArcTail(dec, iterArc);
         spqr_node otherEnd = otherHead == splitNode ? otherTail : otherHead;
         spqr_arc nextArc = getNextNodeArc(dec, iterArc, splitNode);//Need to do this before we modify the arc

         SCIP_Bool isCut = newRow->isArcCut[iterArc];
         SCIP_Bool isMoveColor = newRow->nodeColors[otherEnd] == COLOR_SOURCE;
         SCIP_Bool changeArcEnd = isCut == isMoveColor;
         if( changeArcEnd )
         {
            if( otherHead == splitNode )
            {
               changeArcHead(dec, iterArc, otherHead, newNode);
            } else
            {
               changeArcTail(dec, iterArc, otherTail, newNode);
            }
            if( iterArc == smallMarker )
            {
               smallOtherNode = splitNode;
            }
         }
         newRow->nodeColors[otherEnd] = UNCOLORED;
         //Ugly hack to make sure we can iterate neighbourhood whilst changing arc ends.
         spqr_arc previousArc = iterArc;
         iterArc = nextArc;
         if( iterArc == firstNodeArc )
         {
            break;
         }
         if( changeArcEnd && previousArc == firstNodeArc )
         {
            firstNodeArc = iterArc;
         }
      } while( TRUE ); /*lint !e506*/
   }

   spqr_arc representative = findArcSign(dec, smallMarker).representative;

   newRowInformation->representative = mergeArcSigns(dec, newRowInformation->representative, representative,
                                                     newRow->reducedMembers[smallMember].willBeReversed);

   spqr_node largeMarkerHead = findArcHead(dec, largeMarker);
   spqr_node largeMarkerTail = findArcTail(dec, largeMarker);
   if( findArcSign(dec, largeMarker).reversed )
   {
      spqr_node temp = largeMarkerHead;
      largeMarkerHead = largeMarkerTail;
      largeMarkerTail = temp;
   }
   assert(newRowInformation->head == largeMarkerHead || newRowInformation->head == largeMarkerTail ||
          newRowInformation->tail == largeMarkerHead || newRowInformation->tail == largeMarkerTail);
   spqr_node largeOtherNode = ( newRowInformation->head == largeMarkerHead ||
                                newRowInformation->head == largeMarkerTail )
                              ? newRowInformation->tail
                              : newRowInformation->head;

   spqr_member mergedMember = SPQR_INVALID_MEMBER;
   spqr_node arcNodeOne;
   spqr_node arcNodeTwo;
   spqr_node mergeNodeThree;
   if( largeIsParent )
   {
      SCIP_CALL(
         mergeSplitMemberIntoParent(dec, newRow, smallMemberMember, largeMemberMember, largeMarker, smallMarker, TRUE,
                                    largeOtherNode, smallOtherNode, &mergedMember,
                                    &arcNodeOne,
                                    &arcNodeTwo,
                                    &mergeNodeThree));
   } else
   {
      SCIP_CALL(
         mergeSplitMemberIntoParent(dec, newRow, largeMemberMember, smallMemberMember, smallMarker, largeMarker, TRUE,
                                    smallOtherNode, largeOtherNode, &mergedMember,
                                    &arcNodeOne,
                                    &arcNodeTwo,
                                    &mergeNodeThree));
   }
   newRowInformation->member = mergedMember;

   SCIP_Bool otherIsHead = largeOtherNode == newRowInformation->head;
   SCIP_Bool adjacentToMarkerHead = ( newRowInformation->tail == largeMarkerHead ||
                                      newRowInformation->head == largeMarkerHead );
   if( adjacentToMarkerHead )
   {
      if( otherIsHead )
      {
         newRowInformation->head = mergeNodeThree;
         newRowInformation->tail = arcNodeTwo;
      } else
      {
         newRowInformation->head = arcNodeTwo;
         newRowInformation->tail = mergeNodeThree;
      }
   } else
   {
      if( otherIsHead )
      {
         newRowInformation->head = mergeNodeThree;
         newRowInformation->tail = arcNodeOne;
      } else
      {
         newRowInformation->head = arcNodeOne;
         newRowInformation->tail = mergeNodeThree;
      }
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE splitAndMerge(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id smallMember,
   SCIP_Bool largeIsParent,
   NewRowInformation* const newRowInformation
)
{
   spqr_member member = newRow->reducedMembers[smallMember].member;
   switch( getMemberType(dec, member))
   {
      case SPQR_MEMBERTYPE_RIGID:
      {
         SCIP_CALL(splitAndMergeRigid(dec, newRow, smallMember, largeIsParent, newRowInformation, member));
         break;
      }
      case SPQR_MEMBERTYPE_PARALLEL:
      {
         SCIP_CALL(splitAndMergeParallel(dec, newRow, smallMember, largeIsParent, newRowInformation, member));
         break;
      }
      case SPQR_MEMBERTYPE_SERIES:
      {
         SCIP_CALL(splitAndMergeSeries(dec, newRow, smallMember, largeIsParent, newRowInformation, member));
         break;
      }
      case SPQR_MEMBERTYPE_LOOP:
      case SPQR_MEMBERTYPE_UNASSIGNED:
      default:
         newRow->remainsNetwork = FALSE;
         assert(FALSE);
   }
   return SCIP_OKAY;
}

static SCIP_RETCODE mergeTree(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   reduced_member_id root,
   NewRowInformation* const newRowInformation
)
{
   //We use the same ordering as when finding
   //go to a leaf. We need to start in a leaf to avoid the ambiguity of choosing an orientation
   //in members which have no cut arcs; otherwise, we might choose the wrong one
   reduced_member_id leaf = root;
   while( newRow->reducedMembers[leaf].numChildren != newRow->reducedMembers[leaf].numPropagatedChildren )
   {
      for( int i = 0; i < newRow->reducedMembers[leaf].numChildren; ++i )
      {
         children_idx idx = newRow->reducedMembers[leaf].firstChild + i;
         reduced_member_id child = newRow->childrenStorage[idx];
         if( newRow->reducedMembers[child].type != TYPE_PROPAGATED )
         {
            leaf = child;
            break;
         }
      }
   }
   SCIP_CALL(splitFirstLeaf(dec, newRow, leaf, newRowInformation));

   reduced_member_id baseNode = leaf;
   reduced_member_id nextNode = newRow->reducedMembers[baseNode].parent;

   while( reducedMemberIsValid(nextNode))
   {
      //check this node
      SCIP_CALL(splitAndMerge(dec, newRow, nextNode, FALSE, newRowInformation));

      //Recursively merge the children
      //use a while loop to avoid recursion; we may get stack overflows for large graphs
      MergeTreeCallData * data = newRow->mergeTreeCallData;

      data[0].id = nextNode;
      data[0].currentChild = newRow->reducedMembers[nextNode].firstChild ;
      int depth = 0;
      while(depth >= 0){
         reduced_member_id id = data[depth].id;
         children_idx childidx = data[depth].currentChild;
         if(childidx == newRow->reducedMembers[id].firstChild + newRow->reducedMembers[id].numChildren){
            --depth;
            continue;
         }
         reduced_member_id currentchild = newRow->childrenStorage[childidx];
         data[depth].currentChild += 1;
         //skip this child if we already processed it or it is not merged
         if( currentchild == baseNode || newRow->reducedMembers[currentchild].type == TYPE_PROPAGATED){
            continue;
         }
         SCIP_CALL(splitAndMerge(dec, newRow, currentchild, TRUE, newRowInformation));

         //recursively process the child
         depth += 1;
         assert(depth < newRow->memMergeTreeCallData);
         data[depth].id = currentchild;
         data[depth].currentChild = newRow->reducedMembers[currentchild].firstChild;
      }
      //Move up one layer
      baseNode = nextNode;
      nextNode = newRow->reducedMembers[nextNode].parent;
   }

   return SCIP_OKAY;
}

static SCIP_RETCODE transformComponentRowAddition(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* newRow,
   SPQRRowReducedComponent* component,
   NewRowInformation* const newRowInformation
)
{
   assert(component);
   if( newRow->reducedMembers[component->root].numChildren ==
       newRow->reducedMembers[component->root].numPropagatedChildren )
   {
      //No merging necessary, only a single component
      reduced_member_id reducedMember = component->root;
      assert(reducedMemberIsValid(reducedMember));
      spqr_member member = newRow->reducedMembers[reducedMember].member;
      SPQRMemberType type = getMemberType(dec, member);

      switch( type )
      {
         case SPQR_MEMBERTYPE_RIGID:
            SCIP_CALL(transformSingleRigid(dec, newRow, reducedMember, member, newRowInformation));
            break;
         case SPQR_MEMBERTYPE_PARALLEL:
         {
            SCIP_CALL(transformSingleParallel(dec, newRow, reducedMember, member, newRowInformation));
            break;
         }
         case SPQR_MEMBERTYPE_LOOP:
         case SPQR_MEMBERTYPE_SERIES:
         {

            newRowInformation->member = member;
            cut_arc_id cutArc = newRow->reducedMembers[reducedMember].firstCutArc;
            spqr_arc arc = newRow->cutArcs[cutArc].arc;
            newRowInformation->reversed = arcIsReversedNonRigid(dec, arc) == newRow->cutArcs[cutArc].arcReversed;
            if( type == SPQR_MEMBERTYPE_LOOP )
            {
               if( getNumMemberArcs(dec, member) == 2 )
               {
                  changeLoopToSeries(dec, member);
               } else
               {
               }
            }
            break;
         }
         case SPQR_MEMBERTYPE_UNASSIGNED:
         default:
            assert(FALSE);
            break;
      }

      return SCIP_OKAY;
   }

   SCIP_CALL(mergeTree(dec, newRow, component->root, newRowInformation));

   return SCIP_OKAY;
}


static SCIP_RETCODE SCIPnetrowaddCreate(
   SCIP* scip,
   SCIP_NETROWADD** prowadd
)
{
   assert(scip);
   SCIP_CALL(SCIPallocBlockMemory(scip, prowadd));
   SCIP_NETROWADD* newRow = *prowadd;

   newRow->remainsNetwork = TRUE;

   newRow->reducedMembers = NULL;
   newRow->memReducedMembers = 0;
   newRow->numReducedMembers = 0;

   newRow->reducedComponents = NULL;
   newRow->memReducedComponents = 0;
   newRow->numReducedComponents = 0;

   newRow->memberInformation = NULL;
   newRow->memMemberInformation = 0;
   newRow->numMemberInformation = 0;

   newRow->cutArcs = NULL;
   newRow->memCutArcs = 0;
   newRow->numCutArcs = 0;
   newRow->firstOverallCutArc = INVALID_CUT_ARC;

   newRow->childrenStorage = NULL;
   newRow->memChildrenStorage = 0;
   newRow->numChildrenStorage = 0;

   newRow->newRowIndex = SPQR_INVALID_ROW;

   newRow->newColumnArcs = NULL;
   newRow->newColumnReversed = NULL;
   newRow->memColumnArcs = 0;
   newRow->numColumnArcs = 0;

   newRow->leafMembers = NULL;
   newRow->numLeafMembers = 0;
   newRow->memLeafMembers = 0;

   newRow->decompositionColumnArcs = NULL;
   newRow->decompositionColumnArcReversed = NULL;
   newRow->memDecompositionColumnArcs = 0;
   newRow->numDecompositionColumnArcs = 0;

   newRow->isArcCut = NULL;
   newRow->isArcCutReversed = NULL;
   newRow->memIsArcCut = 0;
   newRow->numIsArcCut = 0;

   newRow->nodeColors = NULL;
   newRow->memNodeColors = 0;

   newRow->articulationNodes = NULL;
   newRow->memArticulationNodes = 0;
   newRow->numArticulationNodes = 0;

   newRow->articulationNodeSearchInfo = NULL;
   newRow->memNodeSearchInfo = 0;

   newRow->crossingPathCount = NULL;
   newRow->memCrossingPathCount = 0;

   newRow->intersectionDFSData = NULL;
   newRow->memIntersectionDFSData = 0;

   newRow->colorDFSData = NULL;
   newRow->memColorDFSData = 0;

   newRow->artDFSData = NULL;
   newRow->memArtDFSData = 0;

   newRow->createReducedMembersCallstack = NULL;
   newRow->memCreateReducedMembersCallstack = 0;

   newRow->intersectionPathDepth = NULL;
   newRow->memIntersectionPathDepth = 0;

   newRow->intersectionPathParent = NULL;
   newRow->memIntersectionPathParent = 0;

   newRow->mergeTreeCallData = NULL;
   newRow->memMergeTreeCallData = 0;

   newRow->temporaryColorArray = NULL;
   newRow->memTemporaryColorArray = 0;

   return SCIP_OKAY;
}

static void SCIPnetrowaddFree(
   SCIP* scip,
   SCIP_NETROWADD** prowadd
)
{
   assert(*prowadd);

   SCIP_NETROWADD* newRow = *prowadd;

   SCIPfreeBlockMemoryArray(scip, &newRow->temporaryColorArray, newRow->memTemporaryColorArray);
   SCIPfreeBlockMemoryArray(scip, &newRow->createReducedMembersCallstack, newRow->memCreateReducedMembersCallstack);
   SCIPfreeBlockMemoryArray(scip, &newRow->artDFSData, newRow->memArtDFSData);
   SCIPfreeBlockMemoryArray(scip, &newRow->colorDFSData, newRow->memColorDFSData);
   SCIPfreeBlockMemoryArray(scip, &newRow->mergeTreeCallData, newRow->memMergeTreeCallData);
   SCIPfreeBlockMemoryArray(scip, &newRow->intersectionDFSData, newRow->memIntersectionDFSData);
   SCIPfreeBlockMemoryArray(scip, &newRow->intersectionPathParent, newRow->memIntersectionPathParent);
   SCIPfreeBlockMemoryArray(scip, &newRow->intersectionPathDepth, newRow->memIntersectionPathDepth);
   SCIPfreeBlockMemoryArray(scip, &newRow->crossingPathCount, newRow->memCrossingPathCount);
   SCIPfreeBlockMemoryArray(scip, &newRow->articulationNodeSearchInfo, newRow->memNodeSearchInfo);
   SCIPfreeBlockMemoryArray(scip, &newRow->articulationNodes, newRow->memArticulationNodes);
   SCIPfreeBlockMemoryArray(scip, &newRow->nodeColors, newRow->memNodeColors);
   SCIPfreeBlockMemoryArray(scip, &newRow->isArcCut, newRow->memIsArcCut);
   SCIPfreeBlockMemoryArray(scip, &newRow->isArcCutReversed, newRow->memIsArcCut);
   SCIPfreeBlockMemoryArray(scip, &newRow->decompositionColumnArcs, newRow->memDecompositionColumnArcs);
   SCIPfreeBlockMemoryArray(scip, &newRow->decompositionColumnArcReversed, newRow->memDecompositionColumnArcs);
   SCIPfreeBlockMemoryArray(scip, &newRow->newColumnArcs, newRow->memColumnArcs);
   SCIPfreeBlockMemoryArray(scip, &newRow->newColumnReversed, newRow->memColumnArcs);
   SCIPfreeBlockMemoryArray(scip, &newRow->childrenStorage, newRow->memChildrenStorage);
   SCIPfreeBlockMemoryArray(scip, &newRow->cutArcs, newRow->memCutArcs);
   SCIPfreeBlockMemoryArray(scip, &newRow->memberInformation, newRow->memMemberInformation);
   SCIPfreeBlockMemoryArray(scip, &newRow->reducedComponents, newRow->memReducedComponents);
   SCIPfreeBlockMemoryArray(scip, &newRow->reducedMembers, newRow->memReducedMembers);
   SCIPfreeBlockMemoryArray(scip, &newRow->leafMembers, newRow->memLeafMembers);
   SCIPfreeBlockMemory(scip, prowadd);
}

static SCIP_RETCODE SCIPnetrowaddCheck(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* rowadd,
   int row,
   const int* nonzcols,
   const double* nonzvals,
   int nnonzs
)
{
   assert(dec);
   assert(rowadd);
   assert(nnonzs == 0 || nonzcols);

   rowadd->remainsNetwork = TRUE;
   cleanUpPreviousIteration(dec, rowadd);

   SCIP_CALL(newRowUpdateRowInformation(dec, rowadd, row, nonzcols, nonzvals, nnonzs));
   SCIP_CALL(constructRowReducedDecomposition(dec, rowadd));
   SCIP_CALL(createReducedDecompositionCutArcs(dec, rowadd));

   SCIP_CALL(determineLeafReducedMembers(dec, rowadd));
   SCIP_CALL(allocateRigidSearchMemory(dec, rowadd));
   SCIP_CALL(allocateTreeSearchMemory(dec, rowadd));
   //Check for each component if the cut arcs propagate through a row tree marker to a cut arc in another component
   //From the leafs inward.
   propagateComponents(dec, rowadd);
   //It can happen that we are not graphic by some of the checked components.
   //In that case, further checking may lead to errors as some invariants that the code assumes will be broken.
   if( rowadd->remainsNetwork )
   {
      for( int i = 0; i < rowadd->numReducedComponents; ++i )
      {
         determineMergeableTypes(dec, rowadd, rowadd->reducedComponents[i].root);
         //exit early if one is not graphic
         if( !rowadd->remainsNetwork )
         {
            break;
         }
      }
   }

   cleanUpRowMemberInformation(rowadd);

   return SCIP_OKAY;
}

static SCIP_RETCODE SCIPnetrowaddAdd(
   SCIP_NETMATDECDATA* dec,
   SCIP_NETROWADD* rowadd
)
{
   assert(rowadd->remainsNetwork);
   if( rowadd->numReducedComponents == 0 )
   {
      spqr_member newMember = SPQR_INVALID_MEMBER;
      SCIP_CALL(createStandaloneParallel(dec, rowadd->newColumnArcs, rowadd->newColumnReversed,
                                         rowadd->numColumnArcs, rowadd->newRowIndex, &newMember));
   } else if( rowadd->numReducedComponents == 1 )
   {
      NewRowInformation information = emptyNewRowInformation();
      SCIP_CALL(transformComponentRowAddition(dec, rowadd, &rowadd->reducedComponents[0], &information));

      if( rowadd->numColumnArcs == 0 )
      {
         spqr_arc rowArc = SPQR_INVALID_ARC;
         SCIP_CALL(createRowArc(dec, information.member, &rowArc, rowadd->newRowIndex, information.reversed));
         if( SPQRnodeIsValid(information.head))
         {
            assert(SPQRnodeIsValid(information.tail));
            assert(SPQRarcIsValid(information.representative));
            setArcHeadAndTail(dec, rowArc, findNode(dec, information.head), findNode(dec, information.tail));
            arcSetRepresentative(dec, rowArc, information.representative);
            arcSetReversed(dec, rowArc, information.reversed != arcIsReversedNonRigid(dec, information.representative));
         }
      } else
      {
         spqr_member new_row_parallel = SPQR_INVALID_MEMBER;
         SCIP_CALL(createConnectedParallel(dec, rowadd->newColumnArcs, rowadd->newColumnReversed, rowadd->numColumnArcs,
                                           rowadd->newRowIndex, &new_row_parallel));
         spqr_arc markerArc = SPQR_INVALID_ARC;
         spqr_arc ignore = SPQR_INVALID_ARC;
         SCIP_CALL(createMarkerPairWithReferences(dec, information.member, new_row_parallel, TRUE,
                                                  information.reversed, FALSE,
                                                  &markerArc, &ignore));
         if( SPQRnodeIsValid(information.head))
         {
            assert(SPQRnodeIsValid(information.tail));
            assert(SPQRarcIsValid(information.representative));
            setArcHeadAndTail(dec, markerArc, findNode(dec, information.head), findNode(dec, information.tail));
            arcSetRepresentative(dec, markerArc, information.representative);
            arcSetReversed(dec, markerArc,
                           information.reversed != arcIsReversedNonRigid(dec, information.representative));
         }
      }
      if( getMemberType(dec, information.member) == SPQR_MEMBERTYPE_LOOP )
      {
         assert(getNumMemberArcs(dec, information.member) == 2 || getNumMemberArcs(dec, information.member) == 3);
         if( getNumMemberArcs(dec, information.member) == 3 )
         {
            changeLoopToSeries(dec, information.member);
         }
      }
   } else
   {
#ifndef NDEBUG
      int numDecComponentsBefore = numConnectedComponents(dec);
#endif
      spqr_member new_row_parallel = SPQR_INVALID_MEMBER;
      SCIP_CALL(createConnectedParallel(dec, rowadd->newColumnArcs, rowadd->newColumnReversed, rowadd->numColumnArcs,
                                        rowadd->newRowIndex, &new_row_parallel));
      for( int i = 0; i < rowadd->numReducedComponents; ++i )
      {
         NewRowInformation information = emptyNewRowInformation();

         SCIP_CALL(transformComponentRowAddition(dec, rowadd, &rowadd->reducedComponents[i], &information));
         if( getMemberType(dec, information.member) == SPQR_MEMBERTYPE_LOOP )
         {
            assert(getNumMemberArcs(dec, information.member) == 1);
            spqr_arc arc = getFirstMemberArc(dec, information.member);
            assert(rowadd->isArcCut[arc]);
            moveArcToNewMember(dec, arc, information.member, new_row_parallel);
            arcSetReversed(dec, arc, rowadd->isArcCutReversed[arc]);
            dec->members[information.member].type = SPQR_MEMBERTYPE_UNASSIGNED;
         } else
         {
            reorderComponent(dec,
                             information.member);//Make sure the new component is the root of the local decomposition tree
            spqr_arc markerArc = SPQR_INVALID_ARC;
            spqr_arc ignore = SPQR_INVALID_ARC;
            SCIP_CALL(createMarkerPairWithReferences(dec, new_row_parallel, information.member, FALSE,
                                                     FALSE, information.reversed,
                                                     &ignore, &markerArc));
            if( SPQRnodeIsValid(information.head))
            {
               assert(SPQRnodeIsValid(information.tail));
               assert(SPQRarcIsValid(information.representative));
               setArcHeadAndTail(dec, markerArc, findNode(dec, information.head), findNode(dec, information.tail));
               arcSetRepresentative(dec, markerArc, information.representative);
               arcSetReversed(dec, markerArc,
                              information.reversed != arcIsReversedNonRigid(dec, information.representative));
            }
         }
      }
      decreaseNumConnectedComponents(dec, rowadd->numReducedComponents - 1);
      assert(numConnectedComponents(dec) == ( numDecComponentsBefore - rowadd->numReducedComponents + 1 ));
   }
   return SCIP_OKAY;
}

static SCIP_Bool SCIPnetrowaddRemainsNetwork(const SCIP_NETROWADD* rowadd)
{
   return rowadd->remainsNetwork;
}

struct SCIP_Netmatdec{
   SCIP_NETMATDECDATA * dec;
   SCIP_NETROWADD * rowadd;
   SCIP_NETCOLADD * coladd;
};

SCIP_RETCODE SCIPnetmatdecCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NETMATDEC**      pdec,               /**< buffer to store pointer to created decomposition */
   int                   nrows,              /**< The maximal number of rows that the decomposition can expect */
   int                   ncols               /**< The maximal number of columns that the decomposition can expect */
)
{
   SCIP_CALL(SCIPallocBlockMemory(scip, pdec));
   SCIP_NETMATDEC* dec = *pdec;
   dec->dec = NULL;
   SCIP_CALL(netMatDecDataCreate(scip, &dec->dec, nrows, ncols));
   dec->rowadd = NULL;
   dec->coladd = NULL;
   return SCIP_OKAY;
}

void SCIPnetmatdecFree(
   SCIP_NETMATDEC**      pdec                /**< pointer to the network matrix decomposition to freed */
)
{
   SCIP_NETMATDEC* dec = *pdec;
   SCIP* scip = dec->dec->env;
   if( dec->coladd != NULL)
   {
      SCIPnetcoladdFree(scip, &dec->coladd);
   }
   if( dec->rowadd != NULL)
   {
      SCIPnetrowaddFree(scip, &dec->rowadd);
   }
   netMatDecDataFree(&dec->dec);
   SCIPfreeBlockMemory(scip,pdec);
}

SCIP_RETCODE SCIPnetmatdecTryAddCol(
   SCIP_NETMATDEC*       dec,                /**< Network matrix decomposition */
   int                   column,             /**< The column to add */
   int*                  nonzrows,           /**< The column's nonzero row indices */
   double*               nonzvals,           /**< The column's nonzero entries */
   int                   nnonzs,             /**< The number of nonzeros in the column */
   SCIP_Bool*            success             /**< Buffer to store whether the column was added */
)
{
   if( dec->coladd == NULL)
   {
      SCIP_CALL(SCIPnetcoladdCreate(dec->dec->env, &dec->coladd));
   }

   SCIP_CALL(SCIPnetcoladdCheck(dec->dec, dec->coladd, column, nonzrows, nonzvals, nnonzs));
   *success = SCIPnetcoladdRemainsNetwork(dec->coladd);
   if( *success )
   {
      SCIP_CALL(SCIPnetcoladdAdd(dec->dec, dec->coladd));
   }
   return SCIP_OKAY;
}

SCIP_RETCODE SCIPnetmatdecTryAddRow(
   SCIP_NETMATDEC*       dec,                /**< Network matrix decomposition */
   int                   row,                /**< The row to add */
   int*                  nonzcols,           /**< The row's nonzero column indices */
   double*               nonzvals,           /**< The row's nonzero entries */
   int                   nnonzs,             /**< The number of nonzeros in the row */
   SCIP_Bool*            success             /**< Buffer to store whether the row was added */
)
{
   if( dec->rowadd == NULL)
   {
      SCIP_CALL(SCIPnetrowaddCreate(dec->dec->env, &dec->rowadd));
   }

   SCIP_CALL(SCIPnetrowaddCheck(dec->dec, dec->rowadd, row, nonzcols, nonzvals, nnonzs));
   *success = SCIPnetrowaddRemainsNetwork(dec->rowadd);
   if( *success )
   {
      SCIP_CALL(SCIPnetrowaddAdd(dec->dec, dec->rowadd));
   }
   return SCIP_OKAY;
}

SCIP_Bool SCIPnetmatdecContainsRow(
   SCIP_NETMATDEC*       dec,                /**< The network matrix decomposition */
   int                   row                 /**< The row index that is checked */
)
{
   return netMatDecDataContainsRow(dec->dec, row);
}

SCIP_Bool SCIPnetmatdecContainsColumn(
   SCIP_NETMATDEC*       dec,                /**< The network matrix decomposition */
   int                   column              /**< The column index that is checked */
)
{
   return netMatDecDataContainsColumn(dec->dec, column);
}

void SCIPnetmatdecRemoveComponent(
   SCIP_NETMATDEC*       dec,                /**< The network matrix decomposition */
   int*                  componentrows,      /**< Pointer to the array of rows to delete */
   int                   nrows,              /**< The number of rows to delete */
   int*                  componentcols,      /**< Pointer to the array of columns to delete */
   int                   ncols               /**< The number of columns to delete */
)
{
   netMatDecDataRemoveComponent(dec->dec,componentrows,nrows,componentcols,ncols);
}

SCIP_Bool SCIPnetmatdecIsMinimal(
   SCIP_NETMATDEC*       dec                 /**< The network matrix decomposition */
)
{
   return netMatDecDataIsMinimal(dec->dec);
}


SCIP_Bool SCIPnetmatdecVerifyCycle(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NETMATDEC*       dec,                /**< The network matrix decomposition */
   int                   column,             /**< The column to check */
   int*                  nonzrowidx,         /**< Array with the column's nonzero row indices */
   double*               nonzvals,           /**< Array with the column's nonzero values */
   int                   nnonzs,             /**< Number of nonzeros in the column */
   int*                  pathrowstorage,     /**< A buffer to hold the computed path's rows. Should have size equal or
                                              *   greater than the number of rows in the decomposition. */
   SCIP_Bool*            pathsignstorage     /**< A buffer to store the computed path's row signs. Should have size
                                              * equal or greater than the number of rows in the decomposition. */
)
{
   return netMatDecDataVerifyCycle(scip,dec->dec,column,nonzrowidx,nonzvals,nnonzs,pathrowstorage,pathsignstorage);
}

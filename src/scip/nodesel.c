/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nodesel.c
 * @brief  datastructures and methods for node selectors
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "nodesel.h"



/** node priority queue data structure */
struct NodePQ
{
   int              len;                /**< number of used element slots */
   int              size;               /**< total number of available element slots */
   void**           slots;              /**< array of element slots */
};

/** node selector */
struct Nodesel
{
   const char*      name;               /**< name of node selector */
   const char*      desc;               /**< description of node selector */
   DECL_NODESELINIT((*nodeselinit));    /**< initialise node selector */
   DECL_NODESELEXIT((*nodeselexit));    /**< deinitialise node selector */
   DECL_NODESELSLCT((*nodeselslct));    /**< node selection method */
   DECL_NODESELCOMP((*nodeselcomp));    /**< node comparison method */
   NODESELDATA*     nodeseldata;        /**< node selector data */
   unsigned int     initialized:1;      /**< is node selector initialized? */
};



/* node priority queue methods */

#define PQ_PARENT(q) (((q)-1)/2)
#define PQ_LEFTCHILD(p) (2*(p)+1)
#define PQ_RIGHTCHILD(p) (2*(p)+2)

static
RETCODE nodepqResize(                   /**< resizes node memory to hold at least the given number of nodes */
   NODEPQ*          nodepq,             /**< pointer to a node priority queue */
   const SET*       set,                /**< global SCIP settings */
   int              minsize             /**< minimal number of storeable nodes */
   )
{
   assert(nodepq != NULL);
   
   if( minsize <= nodepq->size )
      return SCIP_OKAY;

   nodepq->size = SCIPsetCalcTreeGrowSize(set, minsize);
   ALLOC_OKAY( reallocMemoryArray(nodepq->slots, nodepq->size) );

   return SCIP_OKAY;
}

RETCODE SCIPnodepqInit(                 /**< initializes node priority queue */
   NODEPQ**         nodepq              /**< pointer to a node priority queue */
   )
{
   assert(nodepq != NULL);

   ALLOC_OKAY( allocMemory(*nodepq) );
   (*nodepq)->len = 0;
   (*nodepq)->size = 0;
   (*nodepq)->slots = NULL;

   return SCIP_OKAY;
}

void SCIPnodepqFree(                    /**< frees node priority queue, but not the data nodes themselves */
   NODEPQ**         nodepq              /**< pointer to a node priority queue */
   )
{
   assert(nodepq != NULL);

   freeMemoryArrayNull((*nodepq)->slots);
   freeMemory(*nodepq);
}

RETCODE SCIPnodepqInsert(               /**< inserts node into node priority queue */
   NODEPQ*          nodepq,             /**< pointer to a node priority queue */
   const SET*       set,                /**< global SCIP settings */
   NODE*            node                /**< node to be inserted */
   )
{
   SCIP* scip;
   NODESEL* nodesel;
   int pos;
   int cmpresult;

   assert(nodepq != NULL);
   assert(nodepq->len >= 0);
   assert(set != NULL);
   assert(node != NULL);

   scip = set->scip;
   nodesel = set->nodesel;

   CHECK_OKAY( nodepqResize(nodepq, set, nodepq->len+1) );

   /* insert node as leaf in the tree, move it towards the root as long it is better than its parent */
   pos = nodepq->len;
   nodepq->len++;
   while( pos > 0 && (*nodesel->nodeselcomp)(nodesel, scip, node, nodepq->slots[PQ_PARENT(pos)]) < 0 )
   {
      nodepq->slots[pos] = nodepq->slots[PQ_PARENT(pos)];
      pos = PQ_PARENT(pos);
   }
   nodepq->slots[pos] = node;

   return SCIP_OKAY;
}

NODE* SCIPnodepqRemove(                 /**< removes and returns best node from the node priority queue */
   NODEPQ*          nodepq,             /**< pointer to a node priority queue */
   const SET*       set                 /**< global SCIP settings */
   )
{
   SCIP* scip;
   NODESEL* nodesel;
   NODE* root;
   NODE* last;
   int pos;
   int childpos;
   int brotherpos;

   assert(nodepq != NULL);
   assert(nodepq->len >= 0);
   assert(set != NULL);

   if( nodepq->len == 0 )
      return NULL;

   scip = set->scip;
   nodesel = set->nodesel;

   /* remove root node of the tree, move the better child to its parents position until the last node
    * of the queue could be placed in the empty slot */
   root = nodepq->slots[0];
   last = nodepq->slots[nodepq->len-1];
   nodepq->len--;
   pos = 0;
   while( pos < PQ_PARENT(nodepq->len-1) )
   {
      childpos = PQ_LEFTCHILD(pos);
      brotherpos = PQ_RIGHTCHILD(pos);
      if( brotherpos <= nodepq->len
         && (*nodesel->nodeselcomp)(nodesel, scip, nodepq->slots[brotherpos], nodepq->slots[childpos]) < 0 )
         childpos = brotherpos;
      if( (*nodesel->nodeselcomp)(nodesel, scip, last, nodepq->slots[childpos]) <= 0 )
         break;
      nodepq->slots[pos] = nodepq->slots[childpos];
      pos = childpos;
   }
   assert(pos <= nodepq->len);
   nodepq->slots[pos] = last;

   return root;
}

NODE* SCIPnodepqFirst(                  /**< returns the best node of the queue without removing it */
   const NODEPQ*    nodepq              /**< pointer to a node priority queue */
   )
{
   assert(nodepq != NULL);
   assert(nodepq->len >= 0);

   if( nodepq->len == 0 )
      return NULL;

   return nodepq->slots[0];
}

int SCIPnodepqLen(                      /**< returns the number of nodes stored in the node priority queue */
   const NODEPQ*    nodepq              /**< pointer to a node priority queue */
   )
{
   assert(nodepq != NULL);
   assert(nodepq->len >= 0);

   return nodepq->len;
}


/* node selector methods */

RETCODE SCIPnodeselCreate(              /**< creates a node selector */
   NODESEL**        nodesel,            /**< pointer to store node selector */
   const char*      name,               /**< name of node selector */
   const char*      desc,               /**< description of node selector */
   DECL_NODESELINIT((*nodeselinit)),    /**< initialise node selector */
   DECL_NODESELEXIT((*nodeselexit)),    /**< deinitialise node selector */
   DECL_NODESELSLCT((*nodeselslct)),    /**< node selection method */
   DECL_NODESELCOMP((*nodeselcomp)),    /**< node comparison method */
   NODESELDATA*     nodeseldata         /**< node selector data */
   )
{
   ALLOC_OKAY( allocMemory(*nodesel) );
   ALLOC_OKAY( duplicateMemoryArray((*nodesel)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray((*nodesel)->desc, desc, strlen(desc)+1) );
   (*nodesel)->nodeselinit = nodeselinit;
   (*nodesel)->nodeselexit = nodeselexit;
   (*nodesel)->nodeselslct = nodeselslct;
   (*nodesel)->nodeselcomp = nodeselcomp;
   (*nodesel)->nodeseldata = nodeseldata;
   (*nodesel)->initialized = FALSE;

   return SCIP_OKAY;
}
   
RETCODE SCIPnodeselFree(                /**< frees memory of constraint handler */
   NODESEL**        nodesel             /**< pointer to constraint handler data structure */
   )
{
   assert(nodesel != NULL);
   assert(*nodesel != NULL);
   assert(!(*nodesel)->initialized);

   freeMemoryArray((*nodesel)->name);
   freeMemoryArray((*nodesel)->desc);
   freeMemory(*nodesel);

   return SCIP_OKAY;
}

RETCODE SCIPnodeselInit(                /**< initializes node selector */
   NODESEL*         nodesel,            /**< node selector */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(nodesel != NULL);
   assert(scip != NULL);

   if( nodesel->initialized )
   {
      char s[255];
      sprintf(s, "Node selector <%s> already initialized", nodesel->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( (*nodesel->nodeselinit)(nodesel, scip) );
   nodesel->initialized = TRUE;

   return SCIP_OKAY;
}

RETCODE SCIPnodeselExit(                /**< deinitializes node selector */
   NODESEL*         nodesel,            /**< node selector */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(nodesel != NULL);
   assert(scip != NULL);

   if( !nodesel->initialized )
   {
      char s[255];
      sprintf(s, "Node selector <%s> not initialized", nodesel->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   CHECK_OKAY( (*nodesel->nodeselexit)(nodesel, scip) );
   nodesel->initialized = FALSE;

   return SCIP_OKAY;
}

RETCODE SCIPnodeselSelect(              /**< select next node to be processed */
   NODESEL*         nodesel,            /**< node selector */
   SCIP*            scip,               /**< SCIP data structure */   
   NODE**           selnode             /**< pointer to store node to be processed next */
   )
{
   assert(nodesel != NULL);
   assert(scip != NULL);
   assert(selnode != NULL);

   CHECK_OKAY( (*nodesel->nodeselslct)(nodesel, scip, selnode) );

   return SCIP_OKAY;
}

int SCIPnodeselCompare(                 /**< compares two nodes; returns -1/0/+1 if node1 better/equal/worse than node2 */
   NODESEL*         nodesel,            /**< node selector */
   SCIP*            scip,               /**< SCIP data structure */   
   NODE*            node1,              /**< first node to compare */
   NODE*            node2               /**< second node to compare */
   )
{
   assert(nodesel != NULL);
   assert(scip != NULL);
   assert(node1 != NULL);
   assert(node2 != NULL);

   return (*nodesel->nodeselcomp)(nodesel, scip, node1, node2);
}

const char* SCIPnodeselGetName(         /**< gets name of node selector */
   NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return nodesel->name;
}

Bool SCIPnodeselIsInitialized(          /**< is node selector initialized? */
   NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return nodesel->initialized;
}


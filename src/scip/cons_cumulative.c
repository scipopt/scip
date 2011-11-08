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

/**@file   cons_cumulative.c
 * @brief  constraint handler for cumulative constraints
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Jens Schulz
 *
 * Given:
 * - a set of jobs, represented by their integer start time variables \f$S_j\f$, their array of processing times \f$p_j\f$ and of
 *   their demands \f$d_j\f$.
 * - an integer resource capacity \f$C\f$
 *
 * The cumulative constraint ensures that for each point in time \f$t\f$ \f$\sum_{j: S_j \leq t < S_j + p_j} d_j \leq C\f$ holds.
 *
 * Separation:
 * - can be done using binary start time model, see Pritskers, Watters and Wolfe
 * - or by just separating relatively weak cuts on the start time variables
 *
 * Propagation:
 * - time tabling, Klein & Scholl (1999)
 * - Edge-finding from Petr Vilim, adjusted and simplified for dynamic propagation
 *   (2009)
 * - energetic reasoning, see Baptiste, Le Pape, Nuijten (2001)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scipdefplugins.h"
#include "scip/cons_cumulative.h"
#include "scip/cons_linking.h"
#include "scip/cons_knapsack.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "cumulative"
#define CONSHDLR_DESC          "cumulative constraint handler"
#define CONSHDLR_SEPAPRIORITY   2100000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -2040000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -3030000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP

/* default parameter values */
#define DEFAULT_USEBINVARS             FALSE /**< should the binary representation be used? */
#define DEFAULT_LOCALCUTS              FALSE /**< should cuts be added only locally? */
#define DEFAULT_USECOVERCUTS            TRUE /**< should covering cuts be added? */
#define DEFAULT_USECORETIMES            TRUE /**< should core-times be propagated? */
#define DEFAULT_USECORETIMESHOLES      FALSE /**< should core-times be propagated to detect holes? */
#define DEFAULT_USEEDGEFINDING         FALSE /**< should edge finding be used? */
#define DEFAULT_USEENERGETICREASONING  FALSE /**< should energetic reasoning be used? */
#define DEFAULT_CUTSASCONSS             TRUE /**< should the cumulative constraint create the cuts as knapsack constraints? */
#define DEFAULT_MAXNODES               500LL /**< maximum nodes to solve an independent cumulative constraint (-1: unlimited) */


/*
 * Data structures
 */

/** constraint data for cumulative constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< array of variable representing the start time of each job */
   SCIP_CONS**           linkingconss;       /**< array of linking constraints for the integer variables */
   SCIP_ROW**            demandrows;         /**< array of rows of linear relaxation of this problem */
   SCIP_ROW**            scoverrows;         /**< array of rows of small cover cuts of this problem */
   SCIP_ROW**            bcoverrows;         /**< array of rows of big cover cuts of this problem */
   int*                  demands;            /**< array containing corresponding demands */
   int*                  durations;          /**< array containing corresponding durations */
   int                   varssize;           /**< size of the vars-, demands, durations,  and linkingconss-arrays */
   int                   nvars;              /**< number of variables */
   int                   ndemandrows;        /**< number of rows of cumulative constraint for linear relaxation */
   int                   demandrowssize;     /**< size of array rows of demand rows */
   int                   nscoverrows;        /**< number of rows of small cover cuts */
   int                   scoverrowssize;     /**< size of array of small cover cuts */
   int                   nbcoverrows;        /**< number of rows of big cover cuts */
   int                   bcoverrowssize;     /**< size of array of big cover cuts */
   int                   capacity;           /**< available cumulative capacity */
   unsigned int          covercuts:1;        /**< cover cuts are created? */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             usebinvars;         /**< should the binary variables be used? */
   SCIP_Bool             cutsasconss;        /**< should the cumulative constraint create the cuts as knapsack constraints? */
   SCIP_Bool             usecoretimes;       /**< should core-times be propagated? */
   SCIP_Bool             usecoretimesholes;  /**< should core-times be propagated to detect holes? */
   SCIP_Bool             useedgefinding;     /**< should edge finding be used? */
   SCIP_Bool             useenergeticreasoning;/**< should energetic reasoning be used? */
   SCIP_Bool             localcuts;          /**< should cuts be added only locally? */
   SCIP_Bool             usecovercuts;       /**< should covering cuts be added? */

   SCIP_Longint          maxnodes;           /**< maximum nodes spend to solve an independent cumulative constraint (-1: unlimited) */
   SCIP_Longint          lastsepanode;       /**< last node in which separation took place */
};


/*
 * local structure for INFERINFO
 */

/*
 * Propagation rules
 */
enum Proprule
{
   PROPRULE_INVALID              = 0,        /**< propagation was applied without a specific propagation rule */ /*lint !e830*/
   PROPRULE_1_CORETIMES          = 1,        /**< core-time propagator */
   PROPRULE_2_CORETIMEHOLES      = 2,        /**< core-time propagator for holes */
   PROPRULE_3_EDGEFINDING        = 3,        /**< edge-finder */
   PROPRULE_4_ENERGETICREASONING = 4         /**< energetic reasoning */
};
typedef enum Proprule PROPRULE;

/** inference information */
struct InferInfo
{
   union
   {
      struct
      {
         unsigned int    proprule:4;         /**< propagation rule that was applied */
         unsigned int    est:13;             /**< earliest start time of all jobs in conflict set */
         unsigned int    lct:15;             /**< latest completion time of all jobs in conflict set */
      } asbits;
      int                asint;              /**< inference information as a single int value */
   } val;
};
typedef struct InferInfo INFERINFO;

/** converts an integer into an inference information */
static
INFERINFO intToInferInfo(
   int                   i                   /**< integer to convert */
   )
{
   INFERINFO inferinfo;

   inferinfo.val.asint = i;

   return inferinfo;
}

/** converts an inference information into an int */
static
int inferInfoToInt(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return inferinfo.val.asint;
}

/** returns the propagation rule stored in the inference information */
static
PROPRULE inferInfoGetProprule(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return (PROPRULE) inferinfo.val.asbits.proprule;
}

/** returns the earliest start time stored in the inference information */
static
int inferInfoGetEst(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return inferinfo.val.asbits.est;
}

/** returns the latest completion time stored in the inference information */
static
int inferInfoGetLct(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return inferinfo.val.asbits.lct;
}


/** constructs an inference information out of a propagation rule, an earliest start and a latest completion time */
static
INFERINFO getInferInfo(
   PROPRULE              proprule,           /**< propagation rule that deduced the value */
   int                   est,                /**< earliest start time of all jobs in conflict set */
   int                   lct                 /**< latest completion time of all jobs in conflict set */
   )
{
   INFERINFO inferinfo;

   inferinfo.val.asbits.proprule = proprule; /*lint !e641*/
   inferinfo.val.asbits.est = est; /*lint !e732*/
   inferinfo.val.asbits.lct = lct; /*lint !e732*/

   return inferinfo;
}

/*
 * local structure for THETATREE
 */

/** Theta tree node structure */
typedef struct ThetaTreeNode THETATREENODE;
struct ThetaTreeNode
{
   THETATREENODE*        parent;             /**< pointer to the parent node */
   THETATREENODE*        left;               /**< pointer to the left child node */
   THETATREENODE*        right;              /**< pointer to the right child node */
   SCIP_Real             value;              /**< value according to which the tree is ordered */
   SCIP_VAR*             var;                /**< pointer to the variable if node is a leaf or NULL */
   int                   energy;             /**< sum of energies from the leaves in this subtree */
   int                   envelop;            /**< envelop of this subtree */
};

/** Theta tree structure */
struct ThetaTree
{
   THETATREENODE*        superroot;          /**< pointer to the dummy super root node; root is left child */
};
typedef struct ThetaTree THETATREE;

/** returns whether the node is a leaf */
static
SCIP_Bool thetatreeIsLeaf(
   THETATREENODE*        node                /**< node to be evaluated */
   )
{
   assert(node != NULL);
   assert(node->parent != NULL);
   return node->left == NULL && node->right == NULL;
}

/** returns whether the tree is empty */
static
SCIP_Bool thetatreeIsEmpty(
   THETATREE*            tree                /**< tree to be evaluated */
   )
{
   assert(tree != NULL);
   assert(tree->superroot != NULL);

   return tree->superroot->left == NULL;
}

/** returns whether the node is a left child */
static
SCIP_Bool thetatreeIsLeftChild(
   THETATREENODE*        node                /**< node to be evaluated */
   )
{
   assert(node != NULL);
   assert(node->parent != NULL);
   return node->parent->left == node;
}

/** creates an empty theta tree node */
static
SCIP_RETCODE createThetaTreeNode(
   SCIP*                 scip,               /**< SCIP data structure */
   THETATREENODE**       node                /**< node to be created */
   )
{
   SCIP_CALL( SCIPallocMemory(scip, node) );

   (*node)->parent = NULL;
   (*node)->left = NULL;
   (*node)->right = NULL;
   (*node)->value = 0.0;
   (*node)->var = NULL;
   (*node)->energy = 0;
   (*node)->envelop = 0;

   return SCIP_OKAY;
}

/** returns the closest leaf to the given node or NULL if tree is empty */
static
THETATREENODE* findLeafNode(
   THETATREE*            tree,               /**< tree in which the node is searched */
   THETATREENODE*        node                /**< node to be searched for */
   )
{
   THETATREENODE* tmpnode;

   assert(tree != NULL);
   assert(node != NULL);

   if( thetatreeIsEmpty(tree) )
      return NULL;


   tmpnode = tree->superroot->left;

   while( !thetatreeIsLeaf(tmpnode) )
   {
      if( node->value <= tmpnode->value )
         tmpnode = tmpnode->left;
      else
         tmpnode = tmpnode->right;
   }

   return tmpnode;
}

/** updates the envelop and energy on trace */
static
void updateEnvelop(
   THETATREE*            tree,               /**< tree data structure */
   THETATREENODE*        node                /**< node to be updated and its parents */
   )
{
   while( node != tree->superroot )
   {
      assert(node != NULL);
      assert(node->left != NULL);
      assert(node->right != NULL);

      /* update envelop and energy */
      node->envelop = MAX( node->left->envelop + node->right->energy, node->right->envelop);
      node->energy = node->left->energy + node->right->energy;

      /* go to parent */
      node = node->parent;
   }
}


/* inserts the given node into the tree if it is not already inserted */
static
SCIP_RETCODE splitThetaTreeLeaf(
   SCIP*                 scip,               /**< SCIP data structure */
   THETATREE*            tree,               /**< tree data structure */
   THETATREENODE*        splitnode,          /**< node to be split */
   THETATREENODE*        node                /**< node to be inserted */
   )
{
   THETATREENODE* newnode;

   assert(scip != NULL);
   assert(node != NULL);

   /* create a new node as parent of the given ones */
   SCIP_CALL( createThetaTreeNode(scip, &newnode) );
   assert(newnode != NULL);

   newnode->parent = splitnode->parent;

   if( thetatreeIsLeftChild(splitnode) )
   {
      newnode->parent->left = newnode;
   }
   else
   {
      newnode->parent->right = newnode;
   }

   if( node->value < splitnode->value )
   {
      /* node is on the left */
      newnode->left = node;
      newnode->right = splitnode;
      newnode->value = node->value;
   }
   else
   {
      /* split node is on the left */
      newnode->left = splitnode;
      newnode->right = node;
      newnode->value = splitnode->value;
   }

   splitnode->parent = newnode;
   node->parent = newnode;

   updateEnvelop(tree, newnode);

   return SCIP_OKAY;
}

/** creates a theta tree node with variable and sorting value */
static
SCIP_RETCODE thetatreeCreateLeaf(
   SCIP*                 scip,               /**< SCIP data structure */
   THETATREENODE**       node,               /**< node to be created */
   SCIP_VAR*             var,                /**< variable to be stored */
   SCIP_Real             value,              /**< value to be stored */
   int                   energy,             /**< sum of energies from the leaves in this subtree */
   int                   envelop             /**< envelop of this subtree */
   )
{
   assert(var != NULL);

   SCIP_CALL( SCIPallocMemory(scip, node) );

   (*node)->parent = NULL;
   (*node)->left = NULL;
   (*node)->right = NULL;
   (*node)->value = value;
   (*node)->var = var;
   (*node)->energy = energy;
   (*node)->envelop = envelop;

   return SCIP_OKAY;
}


/** creates an empty theta tree */
static
SCIP_RETCODE createThetaTree(
   SCIP*                 scip,               /**< SCIP data structure */
   THETATREE**           tree                /**< tree to be created */
   )
{
   THETATREENODE* node;

   assert(scip != NULL);
   assert(tree != NULL);


   SCIP_CALL( SCIPallocMemory(scip, tree) );

   SCIP_CALL( createThetaTreeNode(scip, &node) );

   (*tree)->superroot = node;

   return SCIP_OKAY;
}

/** frees the theta tree node data structure */
static
SCIP_RETCODE freeThetaTreeNode(
   SCIP*                 scip,               /**< SCIP data structure */
   THETATREENODE**       node                /**< node to be freed */
   )
{
   assert(scip != NULL);
   assert(node != NULL);
   assert(*node != NULL);

   if( (*node)->left != NULL || (*node)->right != NULL )
   {

      if( (*node)->left != NULL )
      {
         SCIP_CALL( freeThetaTreeNode(scip, &((*node)->left) ) );
      }

      if( (*node)->right != NULL )
      {
         SCIP_CALL( freeThetaTreeNode(scip, &((*node)->right) ) );
      }

      (*node)->left = NULL;
      (*node)->right = NULL;
      (*node)->parent = NULL;
      (*node)->var = NULL;
      SCIPfreeMemory(scip, node);
   }

   return SCIP_OKAY;
}

/** frees the theta tree node data structure */
static
SCIP_RETCODE freeThetaTreeLeaf(
   SCIP*                 scip,               /**< SCIP data structure */
   THETATREENODE**       node                /**< node to be freed */
   )
{
   assert(scip != NULL);
   assert(node != NULL);
   assert(*node != NULL);

   assert((*node)->left == NULL);
   assert((*node)->right == NULL);

   (*node)->var = NULL;

   SCIPfreeMemory(scip, node);

   return SCIP_OKAY;
}

/** frees the theta tree data structure */
static
SCIP_RETCODE freeThetaTree(
   SCIP*                 scip,               /**< SCIP data structure */
   THETATREE**           tree                /**< tree to be freed */
   )
{
   assert(scip != NULL);
   assert(tree != NULL);

   if( (*tree)->superroot != NULL )
   {
      SCIP_CALL( freeThetaTreeNode(scip, &((*tree)->superroot) ) );
   }

   SCIPfreeMemory(scip, tree);

   return SCIP_OKAY;
}

/** inserts the given node into the tree if it is not already inserted */
static
SCIP_RETCODE thetatreeInsertLeaf(
   SCIP*                 scip,               /**< SCIP data structure */
   THETATREE*            tree,               /**< tree in which the node is inserted */
   THETATREENODE*        node,               /**< node to be inserted */
   SCIP_Bool*            inserted            /**< pointer to store whether the node could be inserted */
   )
{
   assert(scip != NULL);
   assert(tree != NULL);
   assert(node != NULL);
   assert(inserted != NULL);

   *inserted = FALSE;

   /* if the tree is empty the node will be the root node */
   if( thetatreeIsEmpty(tree) )
   {
      tree->superroot->left = node;
      node->parent = tree->superroot;
      *inserted = TRUE;
      return SCIP_OKAY;
   }
   else
   {
      THETATREENODE* splitleaf;

      /* otherwise find the position to insert the node! */
      splitleaf = findLeafNode(tree, node);

      /* node is already inserted */
      if( node == splitleaf )
         return SCIP_OKAY;

      /* split the 'splitnode' and insert 'node' */
      SCIP_CALL( splitThetaTreeLeaf(scip, tree, splitleaf, node) );
      *inserted = TRUE;
   }

   return SCIP_OKAY;
}

/** return the envelop of the theta tree: \f$max_{\Omega \subseteq \Theta} (C * est_{\Omega} + e_{\Omega})\f$ */
static
int thetaTreeGetEnvelop(
   THETATREE*            tree                /**< tree of which the envelop is returned */
   )
{
   assert(tree != NULL);

   if( thetatreeIsEmpty(tree) )
      return 0;

   return tree->superroot->left->envelop;
}

/*
 * local structure for THETA LAMBDA TREE
 */

/** Theta Lambda tree node structure */
typedef struct TLTreeNode TLTREENODE;
struct TLTreeNode
{
   TLTREENODE*           parent;             /**< pointer to the parent node */
   TLTREENODE*           left;               /**< pointer to the left child node */
   TLTREENODE*           right;              /**< pointer to the right child node */
   SCIP_Real             value;              /**< value according to which the tree is ordered */
   SCIP_VAR*             var;                /**< pointer to the variable if node is a leaf or NULL */
   int                   energy;             /**< sum of energies from the theta-leaves in this subtree */
   int                   envelop;            /**< theta envelop of this subtree */
   int                   energyL;            /**< sum of energies from the lambda-leaves in this subtree */
   int                   envelopL;           /**< lambda envelop of this subtree */
   SCIP_Bool             inTheta;            /**< stores whether this node belongs to the set theta or to lambda */
};

/** Theta lambda tree structure */
struct TLTree
{
   TLTREENODE*           superroot;          /**< pointer to the dummy super root node; root is left child */
};
typedef struct TLTree TLTREE;

/** returns whether the node is a leaf */
static
SCIP_Bool tltreeIsLeaf(
   TLTREENODE*           node                /**< node to be evaluated */
   )
{
   assert(node != NULL);

   return node->left == NULL && node->right == NULL;
}

/** returns whether the node is root node */
static
SCIP_Bool tltreeIsRoot(
   TLTREE*               tree,               /**< tree to be evaluated */
   TLTREENODE*           node                /**< node to be evaluated */
   )
{
   assert(tree != NULL);
   assert(node != NULL);
   assert(tree->superroot != NULL);

   return tree->superroot->left == node;
}

/** returns whether the tree is empty */
static
SCIP_Bool tltreeIsEmpty(
   TLTREE*               tree                /**< tree to be evaluated */
   )
{
   assert(tree != NULL);
   assert(tree->superroot != NULL);

   return tree->superroot->left == NULL;
}

/** returns whether the node is a left child */
static
SCIP_Bool tltreeIsLeftChild(
   TLTREENODE*           node                /**< node to be evaluated */
   )
{
   assert(node != NULL);

   return node->parent->left == node;
}

/** returns whether the node is a right child */
static
SCIP_Bool tltreeIsRightChild(
   TLTREENODE*           node                /**< node to be evaluated */
   )
{
   assert(node != NULL);

   return node->parent->right == node;
}

/** returns the sibling of the node */
static
TLTREENODE* tltreeGetSibling(
   TLTREENODE*           node                /**< node to be evaluated */
   )
{
   assert(node != NULL);
   assert(node->parent != NULL);
   assert(node->parent->left != NULL);
   assert(node->parent->right != NULL);

   if( tltreeIsLeftChild(node) )
      return node->parent->right;

   return node->parent->left;
}

/** creates an empty tltree node */
static
SCIP_RETCODE tltreeCreateNode(
   SCIP*                 scip,               /**< SCIP data structure */
   TLTREENODE**          node                /**< node to be created */
   )
{
   SCIP_CALL( SCIPallocMemory(scip, node) );

   (*node)->parent = NULL;
   (*node)->left = NULL;
   (*node)->right = NULL;
   (*node)->value = 0.0;
   (*node)->var = NULL;
   (*node)->energy = 0;
   (*node)->envelop = 0;
   (*node)->energyL = 0;
   (*node)->envelopL = 0;
   (*node)->inTheta = TRUE;

   return SCIP_OKAY;
}

/** returns the closest leaf to the given node or NULL if tree is empty */
static
TLTREENODE* tltreeFindLeafNode(
   TLTREE*               tree,               /**< tree in which the node is searched */
   TLTREENODE*           node                /**< node to be searched for */
   )
{
   TLTREENODE* tmpnode;

   assert(tree != NULL);
   assert(node != NULL);

   if( tltreeIsEmpty(tree) )
      return NULL;

   tmpnode = tree->superroot->left;

   while( !tltreeIsLeaf(tmpnode) )
   {
      if( node->value <= tmpnode->value )
         tmpnode = tmpnode->left;
      else
         tmpnode = tmpnode->right;
   }

   return tmpnode;
}

/** updates the value of the first parent on the trace which comes from left  */
static
void tltreeUpdateValuesOnTrace(
   TLTREE*               tree,               /**< tree data structure */
   TLTREENODE*           node,               /**< node to be updated or one of its parents */
   SCIP_Real             value               /**< value to be set */
   )
{
   assert(node != NULL);

   while( !tltreeIsRoot(tree, node) )
   {
      if( tltreeIsLeftChild(node) )
      {
         SCIPdebugMessage("update on a trace from %g to %g", node->parent->value, value);
         node->parent->value = value;
         return;
      }
      node = node->parent;
   }
}

/** updates the envelop and energy on trace */
static
void tltreeUpdateEnvelop(
   TLTREE*               tree,               /**< tree data structure */
   TLTREENODE*           node                /**< node to be updated and its parents */
   )
{
   while( node != NULL && node != tree->superroot )
   {
      assert(node != NULL);
      assert(node->left != NULL);
      assert(node->right != NULL);

      /* update envelop and energy */
      node->envelop = MAX( node->left->envelop + node->right->energy, node->right->envelop);
      node->energy = node->left->energy + node->right->energy;

      node->envelopL = MAX( node->left->envelopL + node->right->energy, node->right->envelopL );
      node->envelopL = MAX( node->envelopL , node->left->envelop + node->right->energyL );

      node->energyL = MAX( node->left->energyL + node->right->energy, node->left->energy + node->right->energyL );

      /* negative values are integer min value */
      if( node->envelop < 0 )
         node->envelop = INT_MIN;
      if( node->envelopL < 0 )
         node->envelopL = INT_MIN;
      if( node->energyL < 0 )
         node->energyL = INT_MIN;
      if( node->energy < 0 )
         node->energy = INT_MIN;

      /* go to parent */
      node = node->parent;
   }
}

/** inserts the given node into the tree if it is not already inserted */
static
SCIP_RETCODE tltreeSplitLeaf(
   SCIP*                 scip,               /**< SCIP data structure */
   TLTREE*               tree,               /**< tree data structure */
   TLTREENODE*           splitnode,          /**< node to be split */
   TLTREENODE*           node                /**< node to be inserted */
   )
{
   TLTREENODE* newnode;

   assert(scip != NULL);
   assert(node != NULL);

   /* create a new node as parent of the given ones */
   SCIP_CALL( tltreeCreateNode(scip, &newnode) );
   assert(newnode != NULL);

   newnode->parent = splitnode->parent;

   if( tltreeIsLeftChild(splitnode) )
   {
      newnode->parent->left = newnode;
   }
   else
   {
      newnode->parent->right = newnode;
   }

   if( node->value <= splitnode->value )
   {
      /* node is on the left */
      newnode->left = node;
      newnode->right = splitnode;
      newnode->value = node->value;
   }
   else
   {
      /* split node is on the left */
      newnode->left = splitnode;
      newnode->right = node;
      newnode->value = splitnode->value;
   }

   splitnode->parent = newnode;
   node->parent = newnode;

   tltreeUpdateEnvelop(tree, newnode);

   return SCIP_OKAY;
}

/** creates a theta tree node with variable in theta */
static
SCIP_RETCODE tltreeCreateThetaLeaf(
   SCIP*                 scip,               /**< SCIP data structure */
   TLTREENODE**          node,               /**< node to be created */
   SCIP_VAR*             var,                /**< variable to be stored */
   SCIP_Real             value,              /**< value to be stored */
   int                   energy,             /**< sum of energies from the leaves in this subtree */
   int                   envelop             /**< envelop of this subtree */
   )
{
   assert(var != NULL);

   SCIP_CALL( SCIPallocMemory(scip, node) );

   (*node)->parent = NULL;
   (*node)->left = NULL;
   (*node)->right = NULL;
   (*node)->value = value;
   (*node)->var = var;
   (*node)->energy = energy;
   (*node)->envelop = envelop;
   (*node)->energyL = INT_MIN;
   (*node)->envelopL = INT_MIN;
   (*node)->inTheta = TRUE;

   return SCIP_OKAY;
}

/** creates an empty theta tree */
static
SCIP_RETCODE createTltree(
   SCIP*                 scip,               /**< SCIP data structure */
   TLTREE**              tree                /**< tree to be created */
   )
{
   TLTREENODE* node;

   assert(scip != NULL);
   assert(tree != NULL);


   SCIP_CALL( SCIPallocMemory(scip, tree) );

   SCIP_CALL( tltreeCreateNode(scip, &node) );

   (*tree)->superroot = node;

   return SCIP_OKAY;
}

/** inserts the given node into the tree if it is not already inserted */
static
SCIP_RETCODE tltreeInsertLeaf(
   SCIP*                 scip,               /**< SCIP data structure */
   TLTREE*               tree,               /**< tree in which the node is inserted */
   TLTREENODE*           node,               /**< node to be inserted */
   SCIP_Bool*            inserted            /**< pointer to store whether the node could be inserted */
   )
{
   assert(scip != NULL);
   assert(tree != NULL);
   assert(node != NULL);
   assert(inserted != NULL);

   *inserted = FALSE;

   /* if the tree is empty the node will be the root node */
   if( tltreeIsEmpty(tree) )
   {
      tree->superroot->left = node;
      node->parent = tree->superroot;
      *inserted = TRUE;
      return SCIP_OKAY;
   }
   else
   {
      TLTREENODE* splitleaf;

      /* otherwise find the position to insert the node! */
      splitleaf = tltreeFindLeafNode(tree, node);
      assert(tltreeIsLeaf(splitleaf));
      assert(node != splitleaf);

      /* node is already inserted */
      if( node == splitleaf )
      {
         return SCIP_OKAY;
      }

      /* split the 'splitnode' and insert 'node' */
      SCIP_CALL( tltreeSplitLeaf(scip, tree, splitleaf, node) );
      *inserted = TRUE;
   }

   return SCIP_OKAY;
}

/** creates a full theta lambda tree */
static
SCIP_RETCODE tltreeCreateTree(
   SCIP*                 scip,               /**< SCIP data structure */
   TLTREE**              tree,               /**< pointer to the tree to be created */
   TLTREENODE**          nodes,              /**< leaf nodes to be inserted */
   int*                  perm,               /**< permutation of the nodes to be used */
   int                   nvars               /**< number of leaves */
   )
{
   int j;

   assert(scip != NULL);
   assert(tree != NULL);
   assert(nodes != NULL);
   assert(perm != NULL);

   /* create an empty tree */
   SCIP_CALL( createTltree(scip, tree) );

   for( j = 0; j < nvars; ++j )
   {
      SCIP_Bool inserted;

      SCIP_CALL( tltreeInsertLeaf(scip, *tree, nodes[j], &inserted) );
      assert(inserted);
   }

   return SCIP_OKAY;
}

/** frees the theta lambda tree node data structure, all leaves have to be freed on their own */
static
SCIP_RETCODE freeTltreeNode(
   SCIP*                 scip,               /**< SCIP data structure */
   TLTREENODE**          node                /**< node to be freed */
   )
{
   assert(scip != NULL);
   assert(node != NULL);

   if( tltreeIsLeaf(*node) )
      return SCIP_OKAY;

   if( (*node)->left != NULL || (*node)->right != NULL )
   {
      if( (*node)->left != NULL )
      {
         SCIP_CALL( freeTltreeNode(scip, &((*node)->left) ) );
      }

      if( (*node)->right != NULL )
      {
         SCIP_CALL( freeTltreeNode(scip, &((*node)->right) ) );
      }

      (*node)->left = NULL;
      (*node)->right = NULL;
      (*node)->parent = NULL;
      (*node)->var = NULL;

      SCIPfreeMemory(scip, node);
   }

   return SCIP_OKAY;
}

/** frees the theta lambda tree leaf */
static
SCIP_RETCODE freeTltreeLeaf(
   SCIP*                 scip,               /**< SCIP data structure */
   TLTREENODE**          node                /**< node to be freed */
   )
{
   assert(scip != NULL);
   assert(node != NULL);

   (*node)->left = NULL;
   (*node)->right = NULL;
   (*node)->parent = NULL;
   (*node)->var = NULL;

   SCIPfreeMemory(scip, node);

   return SCIP_OKAY;
}

/** frees the theta tree data structure, BUT: all leaves have to be freed on their own */
static
SCIP_RETCODE freeTltree(
   SCIP*                 scip,               /**< SCIP data structure */
   TLTREE**              tree                /**< tree to be freed */
   )
{
   assert(scip != NULL);
   assert(tree != NULL);

   if( (*tree)->superroot != NULL )
   {
      SCIP_CALL( freeTltreeNode(scip, &((*tree)->superroot) ) );
   }

   SCIPfreeMemory(scip, tree);

   return SCIP_OKAY;
}

/** deletes the given node */
static
SCIP_RETCODE tltreeDeleteLeaf(
   SCIP*                 scip,               /**< SCIP data structure */
   TLTREE*               tree,               /**< tree in which the node is deleted */
   TLTREENODE*           node                /**< node to be deleted */
   )
{

   TLTREENODE* sibling;
   TLTREENODE* parent;
   TLTREENODE* grandparent;

   assert(scip != NULL);
   assert(tree != NULL);
   assert(node != NULL);

   assert(tltreeIsLeaf(node));

   if( tltreeIsRoot(tree, node) )
   {
      node->parent = NULL;
      tree->superroot->left = NULL;
   }

   /* the node belongs to a real subtree */
   sibling = tltreeGetSibling(node);
   assert(sibling != NULL);

   parent = node->parent;
   assert(parent != NULL);

   grandparent = parent->parent;
   assert(grandparent != NULL);

   /* reset parent of sibling */
   sibling->parent = grandparent;

   /* reset child of grandparent to sibling */
   if( tltreeIsLeftChild(parent) )
   {
      grandparent->left = sibling;
   }
   else
   {
      grandparent->right = sibling;

      if( tltreeIsRightChild(parent) )
         tltreeUpdateValuesOnTrace(tree, grandparent, sibling->value);
   }
   tltreeUpdateEnvelop(tree, grandparent);

   SCIPfreeMemory(scip, &parent);

   return SCIP_OKAY;
}

/** return the envelop(theta,lambda) */
static
int tltreeGetEnvelopTL(
   TLTREE*               tree                /**< tree of which the envelop is returned */
   )
{
   assert(tree != NULL);

   if( tltreeIsEmpty(tree) )
      return 0;

   return tree->superroot->left->envelopL;
}

/** transforms the leaf from a theta leaf into a lambda leave */
static
SCIP_RETCODE tltreeTransformLeafTtoL(
   SCIP*                 scip,               /**< SCIP data structure */
   TLTREE*               tree,               /**< tree in which the node is contained as leaf */
   TLTREENODE*           node                /**< node to be transformed */
   )
{
   assert(scip != NULL);
   assert(tree != NULL);
   assert(node != NULL);

   node->envelopL = node->envelop;
   node->energyL = node->energy;

   node->envelop = INT_MIN;
   node->energy = 0;

   node->inTheta = FALSE;

   /* update the energy and envelop values on trace */
   tltreeUpdateEnvelop(tree, node->parent);

   return SCIP_OKAY;
}

/** returns the leaf responsible for the energyL */
static
TLTREENODE* tltreeGetResponsibleLeafEnergyL(
   TLTREENODE*           node                /**< node where the search is continued */
   )
{
   assert(node != NULL);

   if( tltreeIsLeaf(node) )
   {
      assert(!node->inTheta);
      return node;
   }

   if( node->energyL == node->left->energyL + node->right->energy )
      return tltreeGetResponsibleLeafEnergyL(node->left);

   assert(node->energyL == node->left->energy + node->right->energyL);
   return tltreeGetResponsibleLeafEnergyL(node->right);
}

/** returns the leaf responsible for the envelopL */
static
TLTREENODE* tltreeGetResponsibleLeafEnvelopL(
   TLTREENODE*           node                /**< node where the search is continued */
   )
{
   assert(node != NULL);

   if( tltreeIsLeaf(node) )
   {
      assert(!node->inTheta);
      return node;
   }

   if( node->envelopL == node->left->envelopL + node->right->energy )
   {
      return tltreeGetResponsibleLeafEnvelopL(node->left);
   }
   else if( node->envelopL == node->left->envelop + node->right->energyL )
   {
      return tltreeGetResponsibleLeafEnergyL(node->right);
   }

   assert(node->envelopL == node->right->envelopL);

   return tltreeGetResponsibleLeafEnvelopL(node->right);
}

/** returns the leaf responsible for the envelopL */
static
TLTREENODE* tltreeFindResponsibleLeaf(
   TLTREE*               tree                /**< tree to search for responsible leaf */
   )
{
   TLTREENODE* root;

   assert(tree != NULL);

   root = tree->superroot->left;
   assert(root != NULL);

   if( tltreeIsLeaf(root) )
      return NULL;

   return tltreeGetResponsibleLeafEnvelopL(root);
}

/** reports all elements from set theta to generate a conflicting set */
static
void reportSubtreeTheta(
   TLTREENODE*           node,               /**< node whose envelopL needs to be backtraced */
   TLTREENODE***         omegaset,           /**< set to be filled */
   int*                  nelements           /**< pointer to store the number of elements in omegaset */
   )
{
   if( !tltreeIsLeaf(node) )
   {
      reportSubtreeTheta(node->left, omegaset, nelements);
      reportSubtreeTheta(node->right, omegaset, nelements);
   }
   else if( node->inTheta )
   {
      SCIPdebugMessage("add node <%s> as elements %d to omegaset\n", SCIPvarGetName(node->var), *nelements);
      (*omegaset)[*nelements] = node;
      (*nelements)++;
   }
}

/** reports all elements from set theta to generate a conflicting set */
static
void reportEnvelop(
   TLTREENODE*           node,               /**< node whose envelopL needs to be backtraced */
   TLTREENODE***         omegaset,           /**< set to be filled */
   int*                  nelements           /**< pointer to store the number of elements in omegaset */
   )
{
   if( tltreeIsLeaf(node) )
   {
      reportSubtreeTheta(node, omegaset, nelements);
   }
   else if( node->envelop == node->left->envelop + node->right->energy )
   {
      reportEnvelop(node->left, omegaset, nelements);
      reportSubtreeTheta(node->right, omegaset, nelements);
   }
   else
   {
      assert(node->envelop == node->right->envelop);
      reportEnvelop(node->right, omegaset, nelements);
   }
}

/** reports all elements from set theta to generate a conflicting set */
static
void reportEnergyL(
   TLTREENODE*           node,               /**< node whose envelopL needs to be backtraced */
   TLTREENODE***         omegaset,           /**< set to be filled */
   int*                  nelements           /**< pointer to store the number of elements in omegaset */
   )
{
   if( tltreeIsLeaf(node) )
      return;

   if( node->energyL == node->left->energyL + node->right->energy )
   {
      reportEnergyL(node->left, omegaset, nelements);
      reportSubtreeTheta(node->right, omegaset, nelements);
   }
   else
   {
      assert(node->energyL == node->left->energy + node->right->energyL);

      reportSubtreeTheta(node->left, omegaset, nelements);
      reportEnergyL(node->right, omegaset, nelements);
   }
}

/** reports all elements from set theta to generate a conflicting set */
static
void reportEnvelopL(
   TLTREENODE*           node,               /**< node whose envelopL needs to be backtraced */
   TLTREENODE***         omegaset,           /**< set to be filled */
   int*                  nelements           /**< pointer to store the number of elements in omegaset */
   )
{

   /* in a leaf there is no lambda element! */
   if( tltreeIsLeaf(node) )
      return;

   if( node->envelopL == node->left->envelopL + node->right->energy )
   {
      reportEnvelopL(node->left, omegaset, nelements);
      reportSubtreeTheta(node->right, omegaset, nelements);
   }
   else if( node->envelopL == node->left->envelop + node->right->energyL )
   {
      reportEnvelop(node->left, omegaset, nelements);
      reportEnergyL(node->right, omegaset, nelements);
   }
   else
   {
      assert(node->envelopL == node->right->envelopL);

      reportEnvelopL(node->right, omegaset, nelements);
   }
}

/** finds an omega set that leads to a violation
 *  user should take care that this method is only called if the envelop(T,L) > C * lct_j
 *  during edgefinding detection
 *  the array omegaset already needs to be allocated with enough space!
 *  it will be filled with the jobs in non-decreasing order of est_j
 */
static
SCIP_RETCODE tltreeReportOmegaSet(
   SCIP*                 scip,               /**< SCIP data structure */
   TLTREE*               tree,               /**< tree in which the node is contained as leaf */
   TLTREENODE***         omegaset,           /**< set to be filled */
   int*                  nelements           /**< pointer to store the number of elements in omegaset */
   )
{
   assert(scip != NULL);
   assert(tree != NULL);
   assert(omegaset != NULL);
   assert(*omegaset != NULL);
   assert(nelements != NULL);

   *nelements = 0;

   assert(tree->superroot->left->envelopL > 0);

   reportEnvelopL(tree->superroot->left, omegaset, nelements);

   return SCIP_OKAY;
}

/*
 * Local methods
 */

#ifndef NDEBUG
/** converts the given double bound which is integral to an int; in optimized mode the function gets inlined for
 *  performance; in debug mode we check some additional conditions
 */
static
int convertBoundToInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bound               /**< double bound to convert */
   )
{
   assert(SCIPisFeasIntegral(scip, bound));
   assert(SCIPisFeasEQ(scip, bound, (SCIP_Real)(int)(bound + 0.5)));

   return (int)(bound + 0.5);
}
#else
#define convertBoundToInt(x, y) ((int)((y) + 0.5))
#endif

/** creates constraint handler data for cumulative constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   /* create precedence constraint handler data */
   assert(conshdlrdata != NULL);
   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );

   return SCIP_OKAY;
}

/** frees constraint handler data for logic or constraint handler */
static
void conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   SCIPfreeMemory(scip, conshdlrdata);
}

/** prints cumulative constraint to file stream */
static
void consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< cumulative constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int v;

   assert(consdata != NULL);

   /* print coefficients */
   SCIPinfoMessage( scip, file, "cumulative(");

   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(consdata->vars[v] != NULL);
      if( v > 0 )
         SCIPinfoMessage(scip, file, ", ");
      SCIPinfoMessage(scip, file, "<%s>(%d)[%d]", SCIPvarGetName(consdata->vars[v]),
         consdata->durations[v], consdata->demands[v]);
   }
   SCIPinfoMessage(scip, file, ") <= %d", consdata->capacity);
}

/** creates constraint data of cumulative constraint */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to consdata */
   SCIP_VAR**            vars,               /**< array of integer variables */
   SCIP_CONS**           linkingconss,       /**< array of linking constraints for the integer variables, or NULL */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   nvars,              /**< number of variables */
   int                   capacity            /**< available cumulative capacity */
   )
{
   int v;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(vars != NULL || nvars > 0);
   assert(demands != NULL);
   assert(durations != NULL);
   assert(capacity >= 0);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->capacity = capacity;
   (*consdata)->demandrows = NULL;
   (*consdata)->demandrowssize = 0;
   (*consdata)->ndemandrows = 0;
   (*consdata)->scoverrows = NULL;
   (*consdata)->nscoverrows = 0;
   (*consdata)->scoverrowssize = 0;
   (*consdata)->bcoverrows = NULL;
   (*consdata)->nbcoverrows = 0;
   (*consdata)->bcoverrowssize = 0;
   (*consdata)->varssize = nvars;
   (*consdata)->nvars = nvars;
   (*consdata)->covercuts = FALSE;

   if( nvars > 0 )
   {
      assert(vars != NULL); /* for flexelint */

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->demands, demands, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->durations, durations, nvars) );
      (*consdata)->linkingconss = NULL;

#ifndef NDEBUG
      for( v = 0; v < nvars; ++v )
      {
         assert((*consdata)->demands[v] >= 0);
         assert((*consdata)->durations[v] >= 0);
      }
#endif

      if( linkingconss != NULL )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->linkingconss, linkingconss, nvars) );
      }

      /* transform variables, if they are not yet transformed */
      if( SCIPisTransformed(scip) )
      {
         SCIPdebugMessage("get tranformed variables and constraints\n");

         /* get transformed variables and do NOT captures these */
         SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );

         for( v = 0; v < nvars; ++v )
         {
            if( SCIPvarIsActive((*consdata)->vars[v]) && !SCIPdoNotMultaggrVar(scip,  (*consdata)->vars[v]) )
            {
               SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, (*consdata)->vars[v]) );
            }
         }

         if( linkingconss != NULL )
         {
            /* get transformed constraints and captures these */
            SCIP_CALL( SCIPtransformConss(scip, (*consdata)->nvars, (*consdata)->linkingconss, (*consdata)->linkingconss) );

            for( v = 0; v < nvars; ++v )
               assert(SCIPgetConsLinking(scip, (*consdata)->vars[v]) == (*consdata)->linkingconss[v]);
         }
      }
   }
   else
   {
      (*consdata)->vars = NULL;
      (*consdata)->demands = NULL;
      (*consdata)->durations = NULL;
      (*consdata)->linkingconss = NULL;
   }

   return SCIP_OKAY;
}


/** collect linking constraints for each integer variable */
static
SCIP_RETCODE collectLinkingCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< pointer to consdata */
   )
{
   int nvars;
   int v;

   assert(scip != NULL);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   assert(nvars > 0);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->linkingconss, nvars) );

   for( v = 0; v < nvars; ++v )
   {
      SCIP_CONS* cons;
      SCIP_VAR* var;

      var = consdata->vars[v];
      assert(var != NULL);

      SCIPdebugMessage("linking constraint (%d of %d) for variable <%s>\n", v+1, nvars, SCIPvarGetName(var));

      /* create linking constraint if it does not exist yet */
      if( !SCIPexistsConsLinking(scip, var) )
      {
         char name[SCIP_MAXSTRLEN];

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "link(%s)", SCIPvarGetName(var));

         /** creates and captures a linking constraint */
         SCIP_CALL( SCIPcreateConsLinking(scip, &cons, name, var, NULL, 0, 0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE /*TRUE*/, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         consdata->linkingconss[v] = cons;

         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
      else
      {
         consdata->linkingconss[v] = SCIPgetConsLinking(scip, var);
      }

      assert(SCIPexistsConsLinking(scip, var));
      assert(consdata->linkingconss[v] != NULL);
      assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(consdata->linkingconss[v])), "linking") == 0 );
      assert(SCIPgetConsLinking(scip, var) == consdata->linkingconss[v]);
   }

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable in the given cumulative constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_VAR*             var                 /**< variables  */
   )
{
   SCIP_CALL( SCIPunlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

#if 0
/** add rounding locks for the given variable in the given cumulative constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_VAR*             var                 /**< variables  */
   )
{
   SCIP_CALL( SCIPlockVarCons(scip, var, TRUE, TRUE) );

   return SCIP_OKAY;
}
#endif

#ifdef PROFILE_DEBUG
/** output of the given profile */
static
void SCIPprofilePrintOut(
   CUMULATIVEPROFILE*    profile             /**< profile to output */
   )
{
   int t;

   for( t=0; t <profile->ntimepoints; t++ )
   {
      SCIPdebugMessage("tp[%d]: %d -> fc=%d\n", t, profile->timepoints[t], profile-> freecapacities[t]);
   }
}
#endif

/** releases LP rows of constraint data and frees rows array */
static
SCIP_RETCODE consdataFreeRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< constraint data */
   )
{
   int r;

   assert(consdata != NULL);
   assert(*consdata != NULL);

   for( r = 0; r < (*consdata)->ndemandrows; ++r )
   {
      assert((*consdata)->demandrows[r] != NULL);
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->demandrows[r]) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->demandrows, (*consdata)->demandrowssize);

   (*consdata)->ndemandrows = 0;
   (*consdata)->demandrowssize = 0;

   /* free rows of cover cuts */
   for( r = 0; r < (*consdata)->nscoverrows; ++r )
   {
      assert((*consdata)->scoverrows[r] != NULL);
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->scoverrows[r]) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->scoverrows, (*consdata)->scoverrowssize);

   (*consdata)->nscoverrows = 0;
   (*consdata)->scoverrowssize = 0;

   for( r = 0; r < (*consdata)->nbcoverrows; ++r )
   {
      assert((*consdata)->bcoverrows[r] != NULL);
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->bcoverrows[r]) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->bcoverrows, (*consdata)->bcoverrowssize);

   (*consdata)->nbcoverrows = 0;
   (*consdata)->bcoverrowssize = 0;

   (*consdata)->covercuts = FALSE;

   return SCIP_OKAY;
}

/** frees a cumulative constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to linear constraint data */
   )
{
   int varssize;
   int nvars;

   assert(consdata != NULL);
   assert(*consdata != NULL);

   nvars =  (*consdata)->nvars;
   varssize = (*consdata)->varssize;

   if( varssize > 0 )
   {
      int v;

      /* release and free the rows */
      SCIP_CALL( consdataFreeRows(scip, consdata) );

      /* release the linking constraints if they were generated */
      if( (*consdata)->linkingconss != NULL )
      {
         for( v = 0; v < nvars; ++v )
         {
            assert((*consdata)->linkingconss[v] != NULL );
            SCIP_CALL( SCIPreleaseCons(scip, &(*consdata)->linkingconss[v]) );
         }
         SCIPfreeBlockMemoryArray(scip, &(*consdata)->linkingconss, varssize);
      }

      /* free arrays */
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->durations, varssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->demands, varssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, varssize);
   }

   /* free memory */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** check for the given starting time variables with their demands and durations if the cumulative conditions for the
 *  given solution is satisfied
 */
static
SCIP_RETCODE checkCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
   SCIP_Bool*            violated,           /**< pointer to store if the cumulative condition is violated */
   SCIP_CONS*            cons,               /**< constraint which is checked */
   SCIP_Bool             printreason         /**< should the reason for the violation be printed? */
   )
{
   int* startsolvalues;       /* stores when each job is starting */
   int* endsolvalues;         /* stores when each job ends */
   int* startindices;         /* we will sort the startsolvalues, thus we need to know which index of a job it corresponds to */
   int* endindices;           /* we will sort the endsolvalues, thus we need to know which index of a job it corresponds to */

   int freecapacity;
   int curtime;            /* point in time which we are just checking */
   int endindex;           /* index of endsolvalues with: endsolvalues[endindex] > curtime */

   int j;

   assert(scip != NULL);
   assert(violated != NULL);

   (*violated) = FALSE;

   if( nvars == 0 )
      return SCIP_OKAY;

   assert(vars != NULL);
   assert(demands != NULL);
   assert(durations != NULL);

   /* compute time points where we have to check whether capacity constraint is infeasible or not */
   SCIP_CALL( SCIPallocBufferArray(scip, &startsolvalues, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endsolvalues, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startindices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endindices, nvars) );

   /* assign variables, start and endpoints to arrays */
   for( j = 0; j < nvars; ++j )
   {
      /* the constraint of the cumulative constraint handler should be called after the integrality check */
      assert(SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, vars[j])));

      startsolvalues[j] = convertBoundToInt(scip, SCIPgetSolVal(scip, sol, vars[j]));
      startindices[j] = j;

      endsolvalues[j] = startsolvalues[j] + durations[j];
      endindices[j] = j;
   }

   /* sort the arrays not-decreasing according to start solution values and end solution values (and sort the
    * corresponding indices in the same way) */
   SCIPsortIntInt(startsolvalues, startindices, nvars);
   SCIPsortIntInt(endsolvalues, endindices, nvars);

#ifndef NDEBUG
   /* check if the arrays are sorted correctly */
   SCIPdebugMessage("Checking solution <%p> with starting times:\n", (void*)sol);
   SCIPdebugMessage("%i | ", startsolvalues[0]);
   for( j = 1; j < nvars; ++j )
   {
      assert ( startsolvalues[j-1] <= startsolvalues[j] );
      SCIPdebugPrintf("%i | ", startsolvalues[j]);
   }
   SCIPdebugPrintf("\nand end times:\n%i | ",endsolvalues[0]);
   for( j = 1; j < nvars; ++j )
   {
      assert ( endsolvalues[j-1] <= endsolvalues[j] );
      SCIPdebugPrintf("%i | ",endsolvalues[j]);
   }
   SCIPdebugPrintf("\n");
#endif

   endindex = 0;
   freecapacity = capacity;

   /* check each start point of a job whether the capacity is kept or not */
   for( j = 0; j < nvars; ++j )
   {
      curtime = startsolvalues[j];

      /* subtract all capacity needed up to this point */
      freecapacity -= demands[startindices[j]];
      while( j+1 < nvars && startsolvalues[j+1] == curtime )
      {
         j++;
         freecapacity -= demands[startindices[j]];
      }

      /* free all capacity usages of jobs that are no longer running */
      while( endindex < nvars && curtime >= endsolvalues[endindex] )
      {
         freecapacity += demands[endindices[endindex]];
         ++endindex;
      }
      assert(freecapacity <= capacity);

      /* check freecapacity to be smaller than zero */
      if( freecapacity < 0 )
      {
         SCIPdebugMessage("freecapacity = %3d \n", freecapacity);
         (*violated) = TRUE;

         if( printreason )
         {
            int i;

            /* first state the violated constraints */
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );

            /* second state the reason */
            SCIPinfoMessage(scip, NULL,
               "violation: at time point %d available capacity = %d, needed capacity = %d\n",
               curtime, capacity, capacity - freecapacity);

            for(i = 0; i < j; ++i )
            {
               if( startsolvalues[i] + durations[startindices[i]] > curtime )
               {
                  SCIPinfoMessage(scip, NULL, "activity %s, start = %i, duration = %d, demand = %d \n",
                     SCIPvarGetName(vars[startindices[i]]), startsolvalues[i], durations[startindices[i]],
                     demands[startindices[i]]);
               }
            }
         }
         break;
      }
   } /*lint --e{850}*/

   /* free all buffer arrays */
   SCIPfreeBufferArray(scip, &endindices);
   SCIPfreeBufferArray(scip, &startindices);
   SCIPfreeBufferArray(scip, &endsolvalues);
   SCIPfreeBufferArray(scip, &startsolvalues);

   return SCIP_OKAY;
}
/** check if the given constraint is valid; checks each starting point of a job whether the remaining capacity is at
 *  least zero or not. If not (*violated) is set to TRUE
 */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be checked */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   SCIP_Bool*            violated,           /**< pointer to store if the constraint is violated */
   SCIP_Bool             printreason         /**< should the reason for the violation be printed? */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(violated != NULL);

   SCIPdebugMessage("check cumulative constraints <%s>\n", SCIPconsGetName(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check the cumulative condition */
   SCIP_CALL( checkCumulativeCondition(scip, sol, consdata->nvars, consdata->vars,
         consdata->durations, consdata->demands, consdata->capacity, violated, cons, printreason) );

   return SCIP_OKAY;
}

/** checks if the constraint is redundant; that is if its capacity can never be exceeded; therefore we check with
 *  respect to the lower and upper bounds of the integer variables the maximum capacity usage for all event points
 */
static
SCIP_RETCODE consCheckRedundancy(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_Bool*            redundant           /**< pointer to store whether this constraint is redundant */
   )
{

   SCIP_VAR* var;
   int* starttimes;              /* stores when each job is starting */
   int* endtimes;                /* stores when each job ends */
   int* startindices;            /* we will sort the startsolvalues, thus we need to know which index of a job it corresponds to */
   int* endindices;              /* we will sort the endsolvalues, thus we need to know which index of a job it corresponds to */

   int freecapacity;             /* remaining capacity */
   int curtime;                  /* point in time which we are just checking */
   int endindex;                 /* index of endsolvalues with: endsolvalues[endindex] > curtime */

   int j;

   assert(scip != NULL);
   assert(redundant != NULL);

   (*redundant) = TRUE;

   /* if no activities are associated with this cumulative then this constraint is redundant */
   if( nvars == 0 )
      return SCIP_OKAY;

   assert(vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &starttimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endtimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startindices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endindices, nvars) );

   /* assign variables, start and endpoints to arrays */
   for( j = 0; j < nvars; ++j )
   {
      var = vars[j];

      starttimes[j] = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
      startindices[j] = j;

      endtimes[j] =  convertBoundToInt(scip, SCIPvarGetUbLocal(var)) + durations[j];
      endindices[j] = j;
   }

   /* sort the arrays not-decreasing according to startsolvalues and endsolvalues (and sort the indices in the same way) */
   SCIPsortIntInt(starttimes, startindices, nvars);
   SCIPsortIntInt(endtimes, endindices, nvars);

   endindex = 0;
   freecapacity = capacity;

   /* check each start point of a job whether the capacity is violated or not */
   for( j = 0; j < nvars; ++j )
   {
      curtime = starttimes[j];

      /* subtract all capacity needed up to this point */
      freecapacity -= demands[startindices[j]];
      while( j+1 < nvars && starttimes[j+1] == curtime )
      {
         ++j;
         freecapacity -= demands[startindices[j]];
      }

      /* free all capacity usages of jobs the are no longer running */
      while( endtimes[endindex] <= curtime )
      {
         freecapacity += demands[endindices[endindex]];
         ++endindex;
      }
      assert(freecapacity <= capacity);

      /* check freecapacity to be smaller than zero */
      if( freecapacity < 0 )
      {
         (*redundant) = FALSE;
         break;
      }
   } /*lint --e{850}*/

   /* free all buffer arrays */
   SCIPfreeBufferArray(scip, &endindices);
   SCIPfreeBufferArray(scip, &startindices);
   SCIPfreeBufferArray(scip, &endtimes);
   SCIPfreeBufferArray(scip, &starttimes);

   return SCIP_OKAY;
}

/** this method reports all jobs that are running during the given time window (left and right bound) and that exceed
 *  the remaining capacity
 */
static
SCIP_RETCODE analyzeConflictCoreTimesCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_VAR*             infervar,           /**< inference variable */
   int                   leftbound,          /**< left bound of the responsible time window */
   int                   rightbound,         /**< right bound of the responsible time window */
   int                   inferduration,      /**< duration of the inference variable */
   int                   inferdemand,        /**< demand of the inference variable */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_Bool*            success             /**< pointer to store if the conflict was resolved */
   )
{
   SCIP_VAR** corevars;
   SCIP_VAR* var;
   int* startvalues;      /* stores when core of each job is starting */
   int* endvalues;        /* stores when core of each job ends */
   int* startindices;     /* we will sort the startvalues, thus we need to know which index of a job it corresponds to */
   int* endindices;       /* we will sort the endvalues, thus we need to know which index of a job it corresponds to */
   int* coredemands;

   int* conflictids;      /* array where we store job indices of running jobs that are probably in a conflict */
   int nconflictids;

   int j;
   int i;

   int corelb;
   int coreub;
   int freecapacity;       /* remaining capacity */
   int curtime;            /* point in time which we are just checking */
   int endindex;           /* index of endsolvalues with: endsolvalues[endindex] > curtime */

   int ncores;

   assert(nvars > 0);
   assert(leftbound < rightbound);
   assert(success != NULL);
   assert(inferdemand > 0);

   SCIPdebugMessage("analyze reason of '%s' bound change of variable <%s>(%d)[%d], bounds [%d,%d], cap = %d \n",
      boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", SCIPvarGetName(infervar),
      inferduration, inferdemand, leftbound, rightbound, capacity );

   *success = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &corevars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startvalues, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endvalues, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startindices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endindices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coredemands, nvars) );

   ncores = 0;

   /* compute all cores of the variables which lay in the considered time window except the inference variable */
   for( j = 0; j < nvars; ++j )
   {
      var = vars[j];
      assert(var != NULL);

      if( var == infervar )
         continue;

      /* compute cores of jobs; if core overlaps interval of inference variable add this job to the array */
      assert(SCIPisFeasEQ(scip, SCIPvarGetUbAtIndex(var, bdchgidx, TRUE), SCIPvarGetUbAtIndex(var, bdchgidx, FALSE)));
      assert(SCIPisFeasIntegral(scip, SCIPvarGetUbAtIndex(var, bdchgidx, TRUE)));
      assert(SCIPisFeasEQ(scip, SCIPvarGetLbAtIndex(var, bdchgidx, TRUE), SCIPvarGetLbAtIndex(var, bdchgidx, FALSE)));
      assert(SCIPisFeasIntegral(scip, SCIPvarGetLbAtIndex(var, bdchgidx, TRUE)));

      corelb = convertBoundToInt(scip, SCIPvarGetUbAtIndex(var, bdchgidx, TRUE));
      coreub = convertBoundToInt(scip, SCIPvarGetLbAtIndex(var, bdchgidx, TRUE)) + durations[j];

      if( corelb < coreub && leftbound < coreub && rightbound > corelb )
      {
         SCIPdebugMessage("core bounds(%d):%s [%d; %d] <%d>\n",
            j, SCIPvarGetName(var), corelb, coreub, demands[j]);

         corevars[ncores] = var;
         startvalues[ncores] = corelb;
         endvalues[ncores] = coreub;
         coredemands[ncores] = demands[j];
         startindices[ncores] = ncores;
         endindices[ncores] = ncores;
         ncores++;
      }
   }

   /* sort the arrays not-decreasing according to startvalues and endvalues (and sort the indices in the same way) */
   SCIPsortIntInt(startvalues, startindices, ncores);
   SCIPsortIntInt(endvalues, endindices, ncores);

   nconflictids = 0;
   endindex = 0;
   curtime = 0;
   freecapacity = capacity - inferdemand;

   SCIP_CALL( SCIPallocBufferArray(scip, &conflictids, ncores) );

   SCIPdebugMessage("find conflict vars\n");

   /* check each start point of a job whether the capacity is respected  or not */
   for( j = 0; endindex < ncores; ++j ) /*lint !e440*/
   {
      /* subtract all capacity needed up to this point */
      if( j < ncores )
      {
         curtime = startvalues[j];
         freecapacity -= coredemands[startindices[j]];
         conflictids[nconflictids] = startindices[j];
         ++nconflictids;

         SCIPdebugMessage("   start of %d\n", startindices[j]);
         while( j+1 < ncores && startvalues[j+1] <= curtime )
         {
            ++j;
            SCIPdebugMessage("   start of %d\n", startindices[j]);
            freecapacity -= coredemands[startindices[j]];
            conflictids[nconflictids] = startindices[j];
            ++nconflictids;
         }
      }
      else
         curtime = endvalues[endindex];

      SCIPdebugMessage("   endindex=%d, nconflictids=%d\n", endindex, nconflictids);

      /* free all capacity usages of jobs the are no longer running */
      while( endindex < ncores && curtime >= endvalues[endindex] )
      {
         SCIPdebugMessage("   end of %d\n", endindices[endindex]);
         freecapacity += coredemands[endindices[endindex]];

         for( i = 0; i < nconflictids; ++i )
         {
            if( conflictids[i] == endindices[endindex]  )
            {
               conflictids[i] = conflictids[nconflictids-1];
               --nconflictids;
               break;
            }
         }
         ++endindex;
      }
      assert(nconflictids >= 0);

      SCIPdebugMessage("   nconflictids=%d\n", nconflictids);
      SCIPdebugMessage("freecap = %d\n",freecapacity);

      /* check freecapacity to be smaller than zero */
      if( freecapacity < 0 )
      {
         SCIPdebugMessage("freecap = %d\n",freecapacity);

         /* figure out running jobs that are responsible and have to be reported to conflict analysis */
         for( i = 0; i < nconflictids; ++i )
         {
            assert(conflictids[i] < nvars);
            assert(corevars[conflictids[i]] != NULL);
            SCIPdebugMessage("report <%s> with demand %d\n",
               SCIPvarGetName(corevars[conflictids[i]]), coredemands[conflictids[i]]);

            SCIP_CALL( SCIPaddConflictUb( scip, corevars[conflictids[i]], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb( scip, corevars[conflictids[i]], bdchgidx) );

            *success = TRUE;
         }
         nconflictids = 0;
      }
   } /*lint --e{850}*/

   assert(*success);

   /* free all buffer arrays */
   SCIPfreeBufferArray(scip, &conflictids);
   SCIPfreeBufferArray(scip, &coredemands);
   SCIPfreeBufferArray(scip, &endindices);
   SCIPfreeBufferArray(scip, &startindices);
   SCIPfreeBufferArray(scip, &endvalues);
   SCIPfreeBufferArray(scip, &startvalues);
   SCIPfreeBufferArray(scip, &corevars);

   return SCIP_OKAY;
}

/** (only) initialize conflict analysis */
static
SCIP_RETCODE initializeConflictAnalysisCoreTimes(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_VAR*             var,                /**< inference variable */
   int                   leftbound,          /**< left bound of the responsible time window */
   int                   rightbound,         /**< right bound of the responsible time window */
   int                   duration,           /**< duration of the inference variable */
   int                   demand,             /**< demand of the inference variable */
   SCIP_BOUNDTYPE        boundtype          /**< the type of the changed bound (lower or upper bound) */
   )
{
   SCIP_Bool success;

   assert(leftbound < rightbound);

   SCIPdebugMessage("initialize conflict analysis\n");

   /* conflict analysis can only be applied in solving stage and if it is turned on */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPinitConflictAnalysis(scip) );

   /* add lower and upper bound of variable which leads to the infeasibility */
   SCIP_CALL( SCIPaddConflictLb(scip, var, NULL ) );
   SCIP_CALL( SCIPaddConflictUb(scip, var, NULL ) );

   SCIPdebugMessage("add lower and upper bounds of variable <%s>\n", SCIPvarGetName(var));

   SCIP_CALL( analyzeConflictCoreTimesCumulative(scip, nvars, vars, durations, demands, capacity,
         var, leftbound, rightbound, duration, demand, boundtype, NULL, &success) );
   assert(success);

   return SCIP_OKAY;
}

/** this method reports all jobs that are running at time 'timepoint' such that the capacity is exceeded
 *  remaining capacity
 */
static
SCIP_RETCODE analyzeConflictCoreTimesBinvarsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_VAR*             inferbinvar,        /**< inference variable */
   SCIP_VAR*             intvar,             /**< integer variable corresponding to binary variable */
   int                   timepoint,          /**< point in time, where capacity will be exceeded */
   int                   inferdemand,        /**< demand of the inference variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_Bool*            success             /**< pointer to store if the conflict was resolved */
   )
{
   SCIP_VAR** corevars;
   SCIP_VAR* var;
   int* indices;          /* we will sort the demands, thus we need to know which index of a job it corresponds to */
   int* coredemands;

   int j;

   int corelb;
   int coreub;
   int freecapacity;       /* remaining capacity */
   int ncores;

   assert(scip != NULL );
   assert(success != NULL );
   assert(inferdemand > 0 );
   assert(SCIPvarGetType(inferbinvar) == SCIP_VARTYPE_BINARY);
   assert( SCIPvarGetType(intvar) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(intvar) == SCIP_VARTYPE_IMPLINT );

   SCIPdebugMessage("analyze reason of bound change of variable <%s>[%d], cap = %d because of capacity at time %d\n",
      SCIPvarGetName(inferbinvar), inferdemand, capacity, timepoint );

   *success = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &corevars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &indices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coredemands, nvars) );

   ncores = 0;

   /* compute all cores of the variables which lay in the considered time window except the inference variable */
   for( j = 0; j < nvars; ++j )
   {
      var = vars[j];
      assert(var != NULL);

      if( intvar == var )
         continue;

      /* compute cores of jobs; if core overlaps interval of
       * inference variable add this job to the array
       */
      corelb = convertBoundToInt(scip, SCIPvarGetUbAtIndex(var, bdchgidx, TRUE));
      coreub = convertBoundToInt(scip, SCIPvarGetLbAtIndex(var, bdchgidx, TRUE)) + durations[j];

      if( corelb < coreub && timepoint < coreub && timepoint >= corelb )
      {
         SCIPdebugMessage("core bounds(%d):%s [%d; %d] <%d>\n",
            j, SCIPvarGetName(var), corelb, coreub, demands[j]);

         corevars[ncores] = var;
         coredemands[ncores] = demands[j];
         indices[ncores] = ncores;
         ncores++;
      }
   }

   /* sort the arrays not-decreasing according to startvalues and endvalues (and sort the indices in the same way) */
   SCIPsortDownIntInt(coredemands, indices, ncores);

   freecapacity = capacity - inferdemand;

   /* check each start point of a job whether the capacity is respected  or not */
   for( j = 0; (j < ncores) && (freecapacity > 0); ++j )
   {
      freecapacity -= demands[j];

      /* add both bounds of variables with a core at 'timepoint' */
      SCIP_CALL( SCIPaddConflictUb( scip, corevars[indices[j]], bdchgidx) );
      SCIP_CALL( SCIPaddConflictLb( scip, corevars[indices[j]], bdchgidx) );

      *success = TRUE;
   }

   assert(*success);

   /* free all buffer arrays */
   SCIPfreeBufferArray(scip, &coredemands);
   SCIPfreeBufferArray(scip, &indices);
   SCIPfreeBufferArray(scip, &corevars);

   return SCIP_OKAY;
}

/** initialize conflict analysis and analyze conflict */
static
SCIP_RETCODE initializeConflictAnalysisCoreTimesBinvars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_VAR*             binvar,             /**< binary inference variable */
   SCIP_VAR*             intvar,             /**< corresponding starttime variable */
   int                   timepoint,          /**< point in time, where capacity will be exceeded */
   int                   demand              /**< demand of the inference variable */
   )
{
   SCIP_Bool success;

   SCIPdebugMessage("initialize conflict analysis\n");

   /* conflict analysis can only be applied in solving stage and if it is turned on */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPinitConflictAnalysis(scip) );

   /* integer variable is not responsible with its bounds! */
   SCIPdebugMessage("add lower and upper bounds of variable <%s>\n", SCIPvarGetName(binvar));

   SCIP_CALL( analyzeConflictCoreTimesBinvarsCumulative(scip, nvars, vars, durations, demands, capacity,
         binvar, intvar, timepoint, demand, NULL, &success) );
   assert(success);

   return SCIP_OKAY;
}

/** updates the bounds by avoiding core infeasibility */
static
SCIP_RETCODE updateBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_CONS*            cons,               /**< constraint which is propagated */
   CUMULATIVEPROFILE*    profile,            /**< profile to use */
   SCIP_VAR*             var,                /**< the variable the bounds should be updated */
   int                   duration,           /**< the duration of the given variable */
   int                   demand,             /**< the demand of the given variable */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is infeasible */
   )
{
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   INFERINFO inferinfo;

   int lb;
   int ub;
   int newlb;
   int newub;

   lb = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
   ub = convertBoundToInt(scip, SCIPvarGetUbLocal(var));
   assert(lb <= ub);

   /* try to improve the lower bound */
   newlb = SCIPprofileGetEarliestFeasibleStart(profile, lb, ub, duration, demand, &infeasible);
   assert(newlb <= ub || infeasible);

   if( infeasible )
   {
      SCIPdebugMessage("infeasibility detected during change of lower bound of <%s> from %d to %d\n", SCIPvarGetName(var), lb, newlb);

      /* initialize conflict analysis */
      SCIP_CALL( initializeConflictAnalysisCoreTimes(scip, nvars, vars, durations, demands, capacity,
            var, lb, ub+duration, duration, demand, SCIP_BOUNDTYPE_LOWER) );
      *initialized = TRUE;
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   assert(newlb <= ub);
   inferinfo = getInferInfo(PROPRULE_1_CORETIMES, 0, 0);

   SCIP_CALL( SCIPinferVarLbCons(scip, var, (SCIP_Real)newlb, cons, inferInfoToInt(inferinfo), TRUE, &infeasible, &tightened) );
   assert(!infeasible);

   if( tightened )
   {
      SCIPdebugMessage("variable <%s> changes lower bound <%d> -> <%d>\n", SCIPvarGetName(var), lb, newlb);
      ++(*nchgbds);
   }

   /* adjust lower bound */
   lb = MAX(lb,newlb);

   /* get latest start due to cores */
   newub = SCIPprofileGetLatestFeasibleStart(profile, lb, ub, duration, demand, &infeasible);
   assert(newub <= ub);

   /* check whether job fits or not */
   if( infeasible )
   {
      SCIPdebugMessage("infeasibility detected during change of upper bound of <%s> from %d to %d\n", SCIPvarGetName(var), ub, newub);
      /* initialize conflict analysis */
      SCIP_CALL( initializeConflictAnalysisCoreTimes(scip, nvars, vars, durations, demands, capacity,
            var, lb, ub+duration, duration, demand, SCIP_BOUNDTYPE_UPPER) );
      *initialized = TRUE;
      *cutoff = TRUE;
   }
   else
   {
      assert(newub >= lb);
      inferinfo = getInferInfo(PROPRULE_1_CORETIMES, 0, 0);

      /* apply bound change */
      SCIP_CALL( SCIPinferVarUbCons(scip, var, (SCIP_Real)newub, cons, inferInfoToInt(inferinfo), TRUE, &infeasible, &tightened) );
      assert(!infeasible);

      if( tightened )
      {
         SCIPdebugMessage("variable <%s> changes upper bound <%d> -> <%d>\n", SCIPvarGetName(var), ub, newub);
         ++(*nchgbds);
      }
   }

   return SCIP_OKAY;
}

/** a cumulative condition is not satisfied if its capacity is exceeded at a time where jobs cannot be shifted (core)
 *  anymore we build up a cumulative profile of all cores of jobs and try to improve bounds of all jobs
 */
static
SCIP_RETCODE propagateCores(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_CONS*            cons,               /**< constraint which is propagated (needed to SCIPinferVar**Cons()) */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is infeasible */
   )
{
   CUMULATIVEPROFILE* profile;
   SCIP_VAR* var;
   SCIP_Bool* cores;
   SCIP_Bool* fixeds;

   SCIP_Bool infeasible;

   int demand;
   int duration;
   int ncores;
   int j;

   assert(scip != NULL);
   assert(nvars > 0);
   assert(cons != NULL);
   assert(cutoff !=  NULL);
   assert(*cutoff ==  FALSE);
   assert(*initialized ==  FALSE);

   SCIPdebugMessage("check/propagate cores of cumulative condition of constraint <%s>\n", SCIPconsGetName(cons));

   SCIP_CALL(SCIPallocBufferArray(scip, &cores, nvars) );
   SCIP_CALL(SCIPallocBufferArray(scip, &fixeds, nvars) );

   infeasible = FALSE;
   ncores = 0;

   SCIP_CALL( SCIPprofileCreate(scip, &profile, capacity, 4*nvars) );

   /* insert all cores */
   for( j = 0; j < nvars; ++j )
   {
      var = vars[j];
      duration = durations[j];
      demand =  demands[j];
      assert(demand > 0);

      assert(SCIPisFeasIntegral(scip, SCIPvarGetLbLocal(var)));
      assert(SCIPisFeasIntegral(scip, SCIPvarGetUbLocal(var)));

      SCIPprofileInsertCore(scip, profile, var, duration, demand, &cores[j], &fixeds[j], &infeasible);

      if( infeasible )
      {
         SCIPdebugMessage("infeasibility due to cores\n");

         /* initialize conflict analysis */
         SCIP_CALL( initializeConflictAnalysisCoreTimes(scip, nvars, vars, durations, demands, capacity,
               var, convertBoundToInt(scip, SCIPvarGetUbLocal(var)),
               convertBoundToInt(scip, SCIPvarGetLbLocal(var)) + duration, duration, demand,  SCIP_BOUNDTYPE_LOWER) );
         *initialized = TRUE;
         *cutoff = TRUE;
         break;
      }

      if( cores[j] )
         ++ncores;
   }

   if( !(*cutoff) && ncores > 0 )
   {
      /* start checking each job whether bounds can be improved */
      for( j = 0; j < nvars; ++j )
      {
         var = vars[j];
         duration = durations[j];
         demand =  demands[j];
         assert(demand > 0);
         assert(duration > 0);

         if( fixeds[j] )
            continue;

         if( cores[j] )
         {
            SCIPprofileDeleteCore(scip, profile, var, duration, demand, NULL);
         }

         /* try to improve bounds */
         SCIP_CALL( updateBounds(scip, nvars, vars, durations, demands, capacity,
               cons, profile, var, duration, demand, nchgbds, initialized, cutoff) );

         if( *cutoff )
            break;

         /* after updating we might have a new core */
         if( cores[j] || SCIPvarGetLbLocal(var) + duration > convertBoundToInt(scip, SCIPvarGetUbLocal(var)) )
         {
            SCIPprofileInsertCore(scip, profile, var, duration, demand, &cores[j], &fixeds[j], &infeasible);
            assert(cores[j]);
            assert(!infeasible);
         }
      }
   }

   /* free allocated memory */
   SCIPprofileFree(scip, &profile);
   SCIPfreeBufferArray(scip, &fixeds);
   SCIPfreeBufferArray(scip, &cores);

   return SCIP_OKAY;
}

/** updates the binary variables by core-times */
static
SCIP_RETCODE checkForHoles(
   SCIP*                 scip,               /**< SCIP data structure */
   CUMULATIVEPROFILE*    profile,            /**< profile to use */
   SCIP_VAR*             var,                /**< the variable whose binary variables should be updated */
   int                   duration,           /**< the duration of the given variable */
   int                   demand,             /**< the demand of the given variable */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_CONS*            cons,               /**< constraint which is propagated (needed to SCIPinferVar**Cons()) */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is infeasible */
   )
{
   SCIP_VAR** binvars;

   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   int offset;
   int nbinvars;
   int lb;
   int ub;
   int t;
   int pos;

   lb = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
   ub = convertBoundToInt(scip, SCIPvarGetUbLocal(var));
   assert(lb <= ub);

   if( ! SCIPexistsConsLinking(scip, var) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetBinvarsLinking(scip, SCIPgetConsLinking(scip, var), &binvars, &nbinvars) );
   assert(nbinvars > 0 || binvars == NULL);

   if( nbinvars <= 1 )
      return SCIP_OKAY;

   assert(binvars != NULL); /* for flexelint */

   /* check each point in time, whether job can be executed! */
   for( t = lb + 1; t < ub; ++t )
   {
      if( !SCIPprofileIsFeasibleStart(profile, t, duration, demand, &pos) )
      {
         INFERINFO inferinfo;

         offset = SCIPgetOffsetLinking(scip, SCIPgetConsLinking(scip, var));
         assert(binvars[t-offset] != NULL);

         inferinfo = getInferInfo(PROPRULE_2_CORETIMEHOLES, t-offset, profile->timepoints[pos]);
         /* apply bound change */
         SCIP_CALL( SCIPinferVarUbCons(scip, binvars[t-offset], 0.0, cons, inferInfoToInt(inferinfo), TRUE, &infeasible, &tightened) );

         if( infeasible )
         {
            SCIPdebugMessage("infeasibility detected during fixing to zero of var <%s> at time %d not scheduable at %d\n", SCIPvarGetName(binvars[t-offset]), t, profile->timepoints[pos]);

            /* initialize conflict analysis */
            SCIP_CALL( initializeConflictAnalysisCoreTimesBinvars(scip, nvars, vars, durations, demands, capacity,
                  binvars[t-offset], var, profile->timepoints[pos], demand) );
            *initialized = TRUE;
            *cutoff = TRUE;
            return SCIP_OKAY;
         }

         if( tightened )
            ++(*nchgbds);
      }
   }

   return SCIP_OKAY;
}

/** propagates the cores and fixes binary variables, possibly creating holes in the domain */
static
SCIP_RETCODE propagateCoresForHoles(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_CONS*            cons,               /**< constraint which is propagated (needed to SCIPinferVar**Cons()) */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is infeasible */
   )
{
   SCIP_VAR* var;
   CUMULATIVEPROFILE* profile;
   SCIP_Bool* cores;
   SCIP_Bool* fixeds;

   SCIP_Bool infeasible;

   int demand;
   int duration;
   int ncores;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   SCIPdebugMessage("check cores of cumulative constraint <%s>\n", SCIPconsGetName(cons));

   SCIP_CALL(SCIPallocBufferArray(scip, &cores, nvars) );
   SCIP_CALL(SCIPallocBufferArray(scip, &fixeds, nvars) );

   *cutoff =  FALSE;
   infeasible = FALSE;
   ncores = 0;

   SCIP_CALL( SCIPprofileCreate(scip, &profile, capacity, 4*nvars) );

   /* insert all cores */
   for( j = 0; j < nvars; ++j )
   {
      var = vars[j];
      duration = durations[j];
      demand =  demands[j];
      assert(demand > 0);

      assert(SCIPisFeasIntegral(scip, SCIPvarGetLbLocal(var)));
      assert(SCIPisFeasIntegral(scip, SCIPvarGetUbLocal(var)));

      SCIPprofileInsertCore(scip, profile, var, duration, demand, &cores[j], &fixeds[j], &infeasible);
      assert(!infeasible);

      if( cores[j] )
         ++ncores;
   }

   if( !(*cutoff) && ncores > 0 )
   {
      /* start checking each job whether bounds can be improved */
      for( j = 0; j < nvars; ++j )
      {
         var = vars[j];
         duration = durations[j];
         demand =  demands[j];
         assert(demand > 0);
         assert(duration > 0);

         if( fixeds[j] )
            continue;

         if( cores[j] )
         {
            SCIPprofileDeleteCore(scip, profile, var, duration, demand, NULL);
         }

         /* try to improve bounds */
         SCIP_CALL( checkForHoles(scip, profile, var, duration, demand,
               nvars, vars, durations, demands, capacity, cons, nchgbds, initialized, cutoff) );

         if( *cutoff )
            break;

         /* after updating we might have a new core */
         if( cores[j] )
         {
            SCIPprofileInsertCore(scip, profile, var, duration, demand, &cores[j], &fixeds[j], &infeasible);
            assert(cores[j]);
            assert(!infeasible); /* cannot be infeasible; otherwise cutoff in checkForHoles() */
         }
      }
   }

   /* free allocated memory */
   SCIPprofileFree(scip, &profile);
   SCIPfreeBufferArray(scip, &fixeds);
   SCIPfreeBufferArray(scip, &cores);

   return SCIP_OKAY;
}

/** returns TRUE if all demands are smaller than the capacity of the cumulative constraint */
static
SCIP_Bool checkDemands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to be checked */
   )
{
   SCIP_CONSDATA* consdata;
   int capacity;
   int nvars;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   /* if no activities are associated with this cumulative then this constraint is not infeasible, return */
   if( nvars == 0 )
      return TRUE;

   assert(consdata->vars != NULL);
   capacity = consdata->capacity;

   /* check each activity: if demand is larger than capacity the problem is infeasible */
   for( j = 0; j < nvars; ++j )
   {
      if( consdata->demands[j] > capacity )
         return FALSE;
   }

   return TRUE;
}

/** creates covering cuts for jobs violating resource constraints */
static
SCIP_RETCODE createCoverCutsTimepoint(
   SCIP*            scip,                 /**< SCIP data structure */
   SCIP_CONS*       cons,                 /**< constraint to be checked */
   int*             startvalues,          /**< upper bounds on finishing time per job for activities from 0,..., nactivities -1 */
   int              time                  /**< at this point in time covering constraints are valid */
   )
{
   SCIP_VAR** binvars;    /* binary variables of some integer variable */
   SCIP_CONSDATA* consdata;
   SCIP_ROW* row;
   int* flexibleids;
   int* demands;

   char rowname[SCIP_MAXSTRLEN];

   int remainingcap;
   int smallcoversize;    /* size of a small cover */
   int bigcoversize;    /* size of a big cover */
   int nbinvars;
   int offset;
   int nvars;

   int nflexible;
   int D;            /* demand of all jobs up to a certain index */
   int j;
   int i;

   assert(cons != NULL);

   /* get constraint data structure */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL );

   nvars = consdata->nvars;

   /* sort jobs according to demands */
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &flexibleids, nvars) );

   nflexible = 0;
   remainingcap = consdata->capacity;

   /* get all jobs intersecting point 'time' with their bounds */
   for( j = 0; j < nvars; ++j )
   {
      int ub;

      ub = convertBoundToInt(scip, SCIPvarGetUbLocal(consdata->vars[j]));

      /* only add jobs to array if they intersect with point 'time' */
      if( startvalues[j] <= time && ub + consdata->durations[j] > time )
      {
         /* if job is fixed, capacity has to be decreased */
         if( startvalues[j] == ub )
         {
            remainingcap -= consdata->demands[j];
         }
         else
         {
            demands[nflexible] = consdata->demands[j];
            flexibleids[nflexible] = j;
            ++nflexible;
         }
      }
   }
   assert(remainingcap >= 0);

   /* sort demands and job ids */
   SCIPsortIntInt(demands, flexibleids, nflexible);

   /*
    * version 1:
    * D_j := sum_i=0,...,j  d_i, finde j maximal, so dass D_j <= remainingcap
    * erzeuge cover constraint
    *
    */

   /* find maximum number of jobs that can run in parallel (-->coversize = j) */
   D = 0;
   j = 0;

   while( j < nflexible && D <= remainingcap )
   {
      D += demands[j];
      ++j;
   }

   /* j jobs form a conflict, set coversize to 'j - 1' */
   bigcoversize = j-1;
   assert(D > remainingcap);
   assert(bigcoversize < nflexible);

   /* - create a row for all jobs and their binary variables.
    * - at most coversize many binary variables of jobs can be set to one
    */

   /* construct row name */
   (void)SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "capacity_coverbig_%d", time);
   SCIP_CALL( SCIPcreateEmptyRow(scip, &row, rowname, -SCIPinfinity(scip), (SCIP_Real)bigcoversize,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), TRUE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   for( j = 0; j < nflexible; ++j )
   {
      int idx;
      int end;
      int start;
      int lb;
      int ub;

      idx = flexibleids[j];

      /* get and add binvars into var array */
      SCIP_CALL( SCIPgetBinvarsLinking(scip, consdata->linkingconss[idx], &binvars, &nbinvars) );
      assert(nbinvars != 0);
      offset = SCIPgetOffsetLinking(scip, consdata->linkingconss[idx]);

      lb = convertBoundToInt(scip, SCIPvarGetLbLocal(consdata->vars[idx]));
      ub = convertBoundToInt(scip, SCIPvarGetUbLocal(consdata->vars[idx]));
      /* compute start and finishing time */
      start = MAX(lb, time + 1 - consdata->durations[idx]) - offset;
      end =  MIN(time, ub) + 1 - offset;

      /* add all necessary binary variables */
      for( i = start; i < end; ++i )
      {
         assert(i >= 0);
         assert(i < nbinvars);
         assert(binvars[i] != NULL);
         SCIP_CALL( SCIPaddVarToRow(scip, row, binvars[i], 1.0) );
      }
   }

   /* insert and release row */
   SCIP_CALL( SCIPflushRowExtensions(scip, row) );

   if( consdata->bcoverrowssize == 0 )
   {
      consdata->bcoverrowssize = 10;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->bcoverrows, consdata->bcoverrowssize) );
   }
   if( consdata->nbcoverrows == consdata->bcoverrowssize )
   {
      consdata->bcoverrowssize *= 2;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->bcoverrows, consdata->nbcoverrows, consdata->bcoverrowssize) );
   }

   consdata->bcoverrows[consdata->nbcoverrows] = row;
   consdata->nbcoverrows++;

   /*
    * version 2:
    * D_j := sum_i=j,...,0  d_i, finde j minimal, so dass D_j <= remainingcap
    * erzeuge cover constraint und fuege alle jobs i hinzu, mit d_i = d_largest
    */
   /* find maximum number of jobs that can run in parallel (= coversize -1) */
   D = 0;
   j = nflexible -1;
   while( D <= remainingcap )
   {
      assert(j >= 0);
      D += demands[j];
      --j;
   }

   smallcoversize = nflexible - (j + 1) - 1;
   while( j > 0 && demands[j] == demands[nflexible-1] )
      --j;

   assert(smallcoversize < nflexible);

   if( smallcoversize != 1 || smallcoversize != nflexible - (j + 1) - 1 )
   {
      /* construct row name */
      (void)SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "capacity_coversmall_%d", time);
      SCIP_CALL( SCIPcreateEmptyRow(scip, &row, rowname, -SCIPinfinity(scip), (SCIP_Real)smallcoversize,
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), TRUE) );
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

      /* filter binary variables for each unfixed job */
      for( j = j + 1; j < nflexible; ++j )
      {
         int idx;
         int end;
         int start;
         int lb;
         int ub;

         idx = flexibleids[j];

         /* get and add binvars into var array */
         SCIP_CALL( SCIPgetBinvarsLinking(scip, consdata->linkingconss[idx], &binvars, &nbinvars) );
         assert(nbinvars != 0);
         offset = SCIPgetOffsetLinking(scip, consdata->linkingconss[idx]);

         lb = convertBoundToInt(scip, SCIPvarGetLbLocal(consdata->vars[idx]));
         ub = convertBoundToInt(scip, SCIPvarGetUbLocal(consdata->vars[idx]));
         /* compute start and finishing time */
         start = MAX(lb, time + 1 - consdata->durations[idx]) - offset;
         end =  MIN(time, ub) + 1 - offset;

         /* add  all necessary binary variables */
         for( i = start; i < end; ++i )
         {
            assert(i >= 0);
            assert(i < nbinvars);
            assert(binvars[i] != NULL);
            SCIP_CALL( SCIPaddVarToRow(scip, row, binvars[i], 1.0) );
         }
      }

      /* insert and release row */
      SCIP_CALL( SCIPflushRowExtensions(scip, row) );
      if( consdata->scoverrowssize == 0 )
      {
         consdata->scoverrowssize = 10;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->scoverrows, consdata->scoverrowssize) );
      }
      if( consdata->nscoverrows == consdata->scoverrowssize )
      {
         consdata->scoverrowssize *= 2;
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->scoverrows, consdata->nscoverrows, consdata->scoverrowssize) );
      }

      consdata->scoverrows[consdata->nscoverrows] = row;
      consdata->nscoverrows++;
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &flexibleids);
   SCIPfreeBufferArray(scip, &demands);

   return SCIP_OKAY;
}

/** method to construct cover cuts for all points in time */
static
SCIP_RETCODE createCoverCuts(
   SCIP*            scip,                      /**< SCIP data structure */
   SCIP_CONS*       cons                       /**< constraint to be separated */
   )
{
   SCIP_CONSDATA* consdata;

   int* startvalues;        /* stores when each job is starting */
   int* endvalues;          /* stores when each job ends */
   int* startvaluessorted;  /* stores when each job is starting */
   int* endvaluessorted;    /* stores when each job ends */
   int* startindices;     /* we sort the startvalues, so we need to know which index of a job it corresponds to */
   int* endindices;       /* we sort the endvalues, so we need to know which index of a job it corresponds to */

   int nvars;               /* number of jobs for this constraint */
   int freecapacity;        /* remaining capacity */
   int curtime;             /* point in time which we are just checking */
   int endidx;              /* index of endsolvalues with: endsolvalues[endindex] > curtime */

   int j;
   int t;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* if no activities are associated with this resource then this constraint is redundant */
   if( consdata->vars == NULL )
      return SCIP_OKAY;

   nvars = consdata->nvars;

   SCIP_CALL( SCIPallocBufferArray(scip, &startvalues, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endvalues, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startvaluessorted, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endvaluessorted, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startindices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endindices, nvars) );

   /* assign start and endpoints to arrays */
   for( j = 0; j < nvars; ++j )
   {
      startvalues[j] = convertBoundToInt(scip, SCIPvarGetLbLocal(consdata->vars[j]));
      startvaluessorted[j] = startvalues[j];

      endvalues[j] = convertBoundToInt(scip, SCIPvarGetUbLocal(consdata->vars[j])) + consdata->durations[j];
      endvaluessorted[j] = endvalues[j];

      startindices[j] = j;
      endindices[j] = j;
   }

   /* sort the arrays not-decreasing according to startsolvalues and endsolvalues
    * (and sort the indices in the same way) */
   SCIPsortIntInt(startvaluessorted, startindices, nvars);
   SCIPsortIntInt(endvaluessorted, endindices, nvars);

   endidx = 0;
   freecapacity = consdata->capacity;

   /* check each startpoint of a job whether the capacity is kept or not */
   for( j = 0; j < nvars; ++j )
   {
      curtime = startvaluessorted[j];

      /* subtract all capacity needed up to this point */
      freecapacity -= consdata->demands[startindices[j]];

      while( j+1 < nvars && startvaluessorted[j+1] == curtime )
      {
         ++j;
         freecapacity -= consdata->demands[startindices[j]];
      }

      /* free all capacity usages of jobs the are no longer running */
      while( endidx < nvars && curtime >= endvaluessorted[endidx] )
      {
         freecapacity += consdata->demands[endindices[endidx]];
         ++endidx;
      }

      assert(freecapacity <= consdata->capacity);
      assert(endidx <= nvars);

      /* --> endindex - points to the next job which will finish
       *     j        - points to the last job that has been released
       */


      /* check freecapacity to be smaller than zero
       * then we will add cover constraints to the MIP
       */
      if( freecapacity < 0 )
      {
         int nextprofilechange;

         /* we can create covering constraints for each pint in time in interval [curtime; nextprofilechange[ */
         if( j < nvars-1 )
            nextprofilechange = MIN( startvaluessorted[j+1], endvaluessorted[endidx] );
         else
            nextprofilechange = endvaluessorted[endidx];

         for( t = curtime; t < nextprofilechange; ++t )
         {
            SCIPdebugMessage("add cover constraint for time %d\n", curtime);

            /* create covering constraint */
            SCIP_CALL( createCoverCutsTimepoint(scip, cons, startvalues, t)  );

         }
      } /* end if freecapacity > 0 */

   } /*lint --e{850}*/

   consdata->covercuts = TRUE;

   /* free all buffer arrays */
   SCIPfreeBufferArray(scip, &endindices);
   SCIPfreeBufferArray(scip, &startindices);
   SCIPfreeBufferArray(scip, &endvaluessorted);
   SCIPfreeBufferArray(scip, &startvaluessorted);
   SCIPfreeBufferArray(scip, &endvalues);
   SCIPfreeBufferArray(scip, &startvalues);

   return SCIP_OKAY;
}

/** collects all necessary binary variables to represent the jobs which can be active at time point of interest */
static
SCIP_RETCODE collectBinaryVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_VAR***           vars,               /**< pointer to the array to store the binary variables */
   int**                 coefs,              /**< pointer to store the coefficients */
   int*                  nvars,              /**< number if collect binary variables */
   int*                  startindices,       /**< permutation with respect to the start times */
   int                   curtime,            /**< current point in time */
   int                   nstarted,           /**< number of jobs that start before the curtime or at curtime */
   int                   nfinished           /**< number of jobs that finished before curtime or at curtime */
   )
{
   SCIP_VAR** binvars;
   SCIP_VAR* var;
   int nbinvars;
   int nrowvars;
   int startindex;
   int endtime;
   int duration;
   int demand;
   int varidx;
   int offset;
   int minub;
   int size;

   size = 10;
   nrowvars = 0;
   startindex = nstarted - 1;

   SCIP_CALL( SCIPallocBufferArray(scip, vars, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, coefs, size) );

   /* search for the (nstarted - nfinished) jobs which are active at curtime */
   while( nstarted - nfinished > nrowvars )
   {
      /* collect job information */
      varidx = startindices[startindex];
      assert(varidx >= 0 && varidx < consdata->nvars);

      var = consdata->vars[varidx];
      duration = consdata->durations[varidx];
      demand = consdata->demands[varidx];
      assert(var != NULL);

      endtime = convertBoundToInt(scip, SCIPvarGetUbGlobal(var)) + duration;

      /* check the end time of this job is larger than the curtime; in this case the job is still running */
      if( endtime > curtime )
      {
         int tau;  /* counter from curtime - duration + 1 to curtime */

         /* check if the linking constraints exists */
         assert(SCIPexistsConsLinking(scip, var));
         assert(SCIPgetConsLinking(scip, var) != NULL);
         assert(SCIPgetConsLinking(scip, var) == consdata->linkingconss[varidx]);

         /* collect linking constraint information */
         SCIP_CALL( SCIPgetBinvarsLinking(scip, consdata->linkingconss[varidx], &binvars, &nbinvars) );
         offset = SCIPgetOffsetLinking(scip, consdata->linkingconss[varidx]);

         minub = MIN(curtime, endtime - duration);

         for( tau = MAX(curtime - duration + 1, offset); tau <= minub; ++tau )
         {
            assert(tau >= offset && tau < nbinvars + offset);
            assert(binvars[tau-offset] != NULL);

            /* ensure array proper array size */
            if( size == *nvars )
            {
               size *= 2;
               SCIP_CALL( SCIPreallocBufferArray(scip, vars, size) );
               SCIP_CALL( SCIPreallocBufferArray(scip, coefs, size) );
            }

            (*vars)[*nvars] = binvars[tau-offset];
            (*coefs)[*nvars] = demand;
            (*nvars)++;
         }
         nrowvars++;
      }

      startindex--;
   }

   return SCIP_OKAY;
}

/** this method creates a row for time point curtime which insures the capacity restriction of the cumulative constraint */
static
SCIP_RETCODE createCapacityRestriction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be checked */
   int*                  startindices,       /**< permutation with respect to the start times */
   int                   curtime,            /**< current point in time */
   int                   nstarted,           /**< number of jobs that start before the curtime or at curtime */
   int                   nfinished,          /**< number of jobs that finished before curtime or at curtime */
   SCIP_Bool             cutsasconss         /**< should the cumulative constraint create the cuts as constraints? */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** binvars;
   int* coefs;
   int nbinvars;
   char name[SCIP_MAXSTRLEN];
   int capacity;
   int b;

   assert(nstarted > nfinished);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars > 0);

   capacity = consdata->capacity;
   assert(capacity > 0);

   nbinvars = 0;
   SCIP_CALL( collectBinaryVars(scip, consdata, &binvars, &coefs, &nbinvars, startindices, curtime, nstarted, nfinished) );

   /* construct row name */
   (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d[%d]", SCIPconsGetName(cons), nstarted-1, curtime);

   if( cutsasconss )
   {
      SCIP_CONS* lincons;

      /* create linear constraint for the linking between the binary variables and the integer variable */
      SCIP_CALL( SCIPcreateConsKnapsack(scip, &lincons, name, 0, NULL, NULL, (SCIP_Longint)(capacity),
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );

      for( b = 0; b < nbinvars; ++b )
      {
         SCIP_CALL( SCIPaddCoefKnapsack(scip, lincons, binvars[b], (SCIP_Longint)coefs[b]) );
      }

      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   }
   else
   {
      SCIP_ROW* row;

      SCIP_CALL( SCIPcreateEmptyRow(scip, &row, name, -SCIPinfinity(scip), (SCIP_Real)capacity, FALSE, FALSE, SCIPconsIsRemovable(cons)) );
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

      for( b = 0; b < nbinvars; ++b )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, binvars[b], (SCIP_Real)coefs[b]) );
      }

      SCIP_CALL( SCIPflushRowExtensions(scip, row) );
      SCIPdebug( SCIP_CALL(SCIPprintRow(scip, row, NULL)) );

      if( consdata->demandrowssize == 0 )
      {
         consdata->demandrowssize = 10;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->demandrows, consdata->demandrowssize) );
      }
      if( consdata->ndemandrows == consdata->demandrowssize )
      {
         consdata->demandrowssize *= 2;
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->demandrows, consdata->ndemandrows, consdata->demandrowssize) );
      }

      consdata->demandrows[consdata->ndemandrows] = row;
      consdata->ndemandrows++;
   }

   SCIPfreeBufferArrayNull(scip, &binvars);
   SCIPfreeBufferArrayNull(scip, &coefs);

   return SCIP_OKAY;
}

/** initialize the sorted event point arrays */
static
void createSortedEventpoints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int*                  starttimes,         /**< array to store sorted start events */
   int*                  endtimes,           /**< array to store sorted end events */
   int*                  startindices,       /**< permutation with respect to the start times */
   int*                  endindices,         /**< permutation with respect to the end times */
   SCIP_Bool             local               /**< shall local bounds be used */
   )
{
   SCIP_VAR* var;
   int nvars;
   int j;

   nvars = consdata->nvars;

   /* assign variables, start and endpoints to arrays */
   for( j = 0; j < nvars; ++j )
   {
      var = consdata->vars[j];
      if( local )
         starttimes[j] = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
      else
         starttimes[j] = convertBoundToInt(scip, SCIPvarGetLbGlobal(var));

      startindices[j] = j;

      if( local )
         endtimes[j] = convertBoundToInt(scip, SCIPvarGetUbLocal(var)) + consdata->durations[j];
      else
         endtimes[j] = convertBoundToInt(scip, SCIPvarGetUbGlobal(var)) + consdata->durations[j];

      endindices[j] = j;
   }

   /* sort the arrays not-decreasing according to startsolvalues and endsolvalues (and sort the indices in the same way) */
   SCIPsortIntInt(starttimes, startindices, nvars);
   SCIPsortIntInt(endtimes, endindices, nvars);
}

/** remove the capacity requirments for all job which start at the curtime */
static
void subtractStartingJobDemands(
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   curtime,            /**< current point in time */
   int*                  starttimes,         /**< array of start times */
   int*                  startindices,       /**< permutation with respect to the start times */
   int*                  freecapacity,       /**< pointer to store the resulting free capacity */
   int*                  idx,                /**< pointer to index in start time array */
   int                   nvars               /**< number of vars in array of starttimes and startindices */
   )
{
   assert(idx != NULL);
#ifdef SCIP_DEBUG
   int oldidx;
   oldidx = *idx;
#endif
   assert(starttimes != NULL);
   assert(starttimes != NULL);
   assert(freecapacity != NULL);
   assert(starttimes[*idx] == curtime);
   assert(consdata->demands != NULL);
   assert(freecapacity != idx);

   /* subtract all capacity needed up to this point */
   (*freecapacity) -= consdata->demands[startindices[*idx]];

   while( (*idx)+1 < nvars && starttimes[(*idx)+1] == curtime )
   {
      ++(*idx);
      (*freecapacity) -= consdata->demands[startindices[(*idx)]];
      assert(freecapacity != idx);
   }
#ifdef SCIP_DEBUG
   assert(oldidx <= *idx);
#endif
}

/** add the capacity requirments for all job which end at the curtime */
static
void addEndingJobDemands(
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   curtime,            /**< current point in time */
   int*                  endtimes,           /**< array of end times */
   int*                  endindices,         /**< permutation with respect to the end times */
   int*                  freecapacity,       /**< pointer to store the resulting free capacity */
   int*                  idx,                /**< pointer to index in end time array */
   int                   nvars               /**< number of vars in array of starttimes and startindices */
   )
{
#ifdef SCIP_DEBUG
   int oldidx;
   oldidx = *idx;
#endif

   /* free all capacity usages of jobs the are no longer running */
   while( endtimes[*idx] <= curtime && *idx < nvars)
   {
      (*freecapacity) += consdata->demands[endindices[*idx]];
      ++(*idx);
   }

#ifdef SCIP_DEBUG
   assert(oldidx <= *idx);
#endif
}

/** this method checks how many cumulatives can run at most at one time if this is greater than the capacity it creates
 *  row
 */
static
SCIP_RETCODE consCapacityConstraintsFinder(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be checked */
   SCIP_Bool             cutsasconss         /**< should the cumulative constraint create the cuts as constraints? */
   )
{
   SCIP_CONSDATA* consdata;

   int* starttimes;         /* stores when each job is starting */
   int* endtimes;           /* stores when each job ends */
   int* startindices;       /* we will sort the startsolvalues, thus we need to know which index of a job it corresponds to */
   int* endindices;         /* we will sort the endsolvalues, thus we need to know which index of a job it corresponds to */

   int nvars;               /* number of activities for this constraint */
   int freecapacity;        /* remaining capacity */
   int curtime;             /* point in time which we are just checking */
   int endindex;            /* index of endsolvalues with: endsolvalues[endindex] > curtime */

   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;


   /* if no activities are associated with this cumulative then this constraint is redundant */
   if( nvars == 0 )
      return SCIP_OKAY;

   assert(consdata->vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &starttimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endtimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startindices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endindices, nvars) );

   SCIPdebugMessage("create sorted event points for cumulative constraint <%s> with %d jobs\n",
      SCIPconsGetName(cons), nvars);

   /* create event point arrays */
   createSortedEventpoints(scip, consdata,starttimes, endtimes, startindices, endindices, FALSE);

   endindex = 0;
   freecapacity = consdata->capacity;

   /* check each startpoint of a job whether the capacity is kept or not */
   for( j = 0; j < nvars; ++j )
   {
      curtime = starttimes[j];
      SCIPdebugMessage("look at %d-th job with start %d\n", j, curtime);

      /* remove the capacity requirments for all job which start at the curtime */
      subtractStartingJobDemands(consdata, curtime, starttimes, startindices, &freecapacity, &j, nvars);

      /* add the capacity requirments for all job which end at the curtime */
      addEndingJobDemands(consdata, curtime, endtimes, endindices, &freecapacity, &endindex, nvars);

      assert(freecapacity <= consdata->capacity);
      assert(endindex <= nvars);

      /* endindex - points to the next job which will finish */
      /* j - points to the last job that has been released */

      /* if free capacity is smaller than zero, then add rows to the LP */
      if( freecapacity < 0 )
      {
         int nextstarttime;
         int t;

         /* step forward until next job is released and see whether capacity constraint is met or not */
         if( j < nvars-1 )
            nextstarttime = starttimes[j+1];
         else
            nextstarttime = endtimes[nvars-1];

         /* create capacity restriction row for current event point */
         SCIP_CALL( createCapacityRestriction(scip, cons, startindices, curtime, j+1, endindex, cutsasconss) );

         /* create for all points in time between the current event point and next start event point a row if the free
          * capacity is still smaller than zero  */
         for( t = curtime+1 ; t < nextstarttime; ++t )
         {
            /* add the capacity requirments for all job which end at the curtime */
            addEndingJobDemands(consdata, t, endtimes, endindices, &freecapacity, &endindex, nvars);

            if( freecapacity < 0 )
            {
               /* add constraint */
               SCIPdebugMessage("add capacity constraint at time %d\n", t);

               /* create capacity restriction row */
               SCIP_CALL( createCapacityRestriction(scip, cons, startindices, t, j+1, endindex, cutsasconss) );
            }
            else
               break;
         }
      }
   } /*lint --e{850}*/

   /* free all buffer arrays */
   SCIPfreeBufferArray(scip, &endindices);
   SCIPfreeBufferArray(scip, &startindices);
   SCIPfreeBufferArray(scip, &endtimes);
   SCIPfreeBufferArray(scip, &starttimes);

   return SCIP_OKAY;
}

/** creates LP rows corresponding to cumulative constraint; therefore, check each point in time if the maximal needed
 *  capacity is larger than the capacity of the cumulative constraint
 *  - for each necessary point in time:
 *
 *    sum_j sum_t demand_j * x_{j,t} <= capacity
 *
 *    where x(j,t) is the binary variables of job j at time t
 */
static
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   SCIP_Bool             cutsasconss         /**< should the cumulative constraint create the cuts as constraints? */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->demandrows == NULL);
   assert(consdata->ndemandrows == 0);

   /* collect the linking constraints */
   if( consdata->linkingconss == NULL )
   {
      SCIP_CALL( collectLinkingCons(scip, consdata) );
   }

   SCIP_CALL( consCapacityConstraintsFinder(scip, cons, cutsasconss) );

   /* switch of separation for the cumulative constraint if linear constraints are add as cuts */
   if( cutsasconss )
   {
      if( SCIPconsIsInitial(cons) )
      {
         SCIP_CALL( SCIPsetConsInitial(scip, cons, FALSE) );
      }
      if( SCIPconsIsSeparated(cons) )
      {
         SCIP_CALL( SCIPsetConsSeparated(scip, cons, FALSE) );
      }
      if( SCIPconsIsEnforced(cons) )
      {
         SCIP_CALL( SCIPsetConsEnforced(scip, cons, FALSE) );
      }
   }

   return SCIP_OKAY;
}

/** adds linear relaxation of cumulative constraint to the LP */
static
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   SCIP_Bool             cutsasconss         /**< should the cumulative constraint create the cuts as constraints? */
   )
{
   SCIP_CONSDATA* consdata;
   int r;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->demandrows == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons, cutsasconss) );
   }

   for( r = 0; r < consdata->ndemandrows; ++r )
   {
      if( !SCIProwIsInLP(consdata->demandrows[r]) )
      {
         assert(consdata->demandrows[r] != NULL);
         SCIP_CALL( SCIPaddCut(scip, NULL, consdata->demandrows[r], FALSE) );
      }
   }

   return SCIP_OKAY;
}

/** repropagation of energetic reasoning algorithm */
static
SCIP_RETCODE analyzeConflictEnergeticReasoning(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   SCIP_VAR*             infervar,           /**< variable whose bound change is to be explained */
   INFERINFO             inferinfo,          /**< inference info containing position of correct bdchgids */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_Bool*            success             /**< pointer to store if we could explain the bound change */
   )
{
   int est;
   int lct;

   int j;

   assert(scip != NULL);
   assert(nvars > 0);
   assert(inferInfoGetProprule(inferinfo) == PROPRULE_4_ENERGETICREASONING);

   SCIPdebugMessage("repropagate energetic reasoning for variable <%s>\n",
      infervar == NULL ? "null" : SCIPvarGetName(infervar));

   *success = FALSE;

   est = inferInfoGetEst(inferinfo);
   lct = inferInfoGetLct(inferinfo);
   assert(est < lct);

   /* collect the current lower bound of the start variables */
   for( j = 0; j < nvars; ++j )
   {
      SCIP_VAR* var;
      var = vars[j];

      if( var == infervar )
         continue;

      /* report all jobs with non-empty intersection with [est_omega; lct_omega]*/
      if( convertBoundToInt(scip, SCIPvarGetLbAtIndex(var, bdchgidx, FALSE)) + durations[j] >= est
         && convertBoundToInt(scip, SCIPvarGetUbAtIndex(var, bdchgidx, FALSE)) <= lct )
      {
         SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx) );
         *success = TRUE;
      }
   }

   if( !(*success) )
   {
      SCIPinfoMessage(scip, NULL, "could not resolve conflict from energetic reasoning\n");
      SCIPABORT();
   }

   return SCIP_OKAY;
}

/** initialize conflict analysis and analyze conflict */
static
SCIP_RETCODE initializeConflictAnalysisEnergeticReasoning(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   SCIP_VAR*             infervar,           /**< inference variable */
   INFERINFO             inferinfo           /**< inference information */
   )
{
   SCIP_Bool success;

   SCIPdebugMessage("initialize conflict analysis for energetic reasoning\n");

   /* conflict analysis can only be applied in solving stage and if it is turned on */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   assert(inferInfoGetProprule(inferinfo) == PROPRULE_4_ENERGETICREASONING);

   SCIP_CALL( SCIPinitConflictAnalysis(scip) );

   /* add lower and upper bound of variable which lead to the infeasibility */
   if( infervar != NULL )
   {
      SCIP_CALL( SCIPaddConflictLb(scip, infervar, NULL ) );
      SCIP_CALL( SCIPaddConflictUb(scip, infervar, NULL ) );

      SCIPdebugMessage("add lower and upper bounds of variable <%s>\n", SCIPvarGetName(infervar));
   }

   SCIP_CALL( analyzeConflictEnergeticReasoning(scip, nvars, vars, durations,
         infervar, inferinfo, NULL, &success) );
   assert(success);

   return SCIP_OKAY;
}

/** computes the energy in the interval [est,lct] of the given variable if the corresponding job is right shifted */
static
int getVarRightEnergy(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   SCIP_HASHMAP*         varhashmap,         /**< hashmap from variables to their indices */
   SCIP_VAR*             var,                /**< variable whose energy is reported */
   int                   est,                /**< left bound of interval */
   int                   lct                 /**< right bound of interval */
   )
{
   int energy;
   int lst_j;
   int min;
   int j;

   assert(scip != NULL);
   assert(varhashmap != NULL);
   assert(var != NULL);
   assert(est < lct);

   SCIPdebugMessage("perform energetic reasoning\n");

   j = (int)(size_t)SCIPhashmapGetImage(varhashmap, var);

   /* compute variable's bounds */
   lst_j = convertBoundToInt(scip, SCIPvarGetUbLocal(var));

   min = MIN(lct - est, lct - lst_j);

   energy =  MAX(0, MIN(min, durations[j])) * demands[j];

   return energy;
}

/** computes the energy in the interval [est,lct] of the given variable if the corresponding job is left shifted */
static
int getVarLeftEnergy(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   SCIP_HASHMAP*         varhashmap,         /**< hashmap from variables to their indices */
   SCIP_VAR*             var,                /**< variable whose energy is reported */
   int                   est,                /**< left bound of interval */
   int                   lct                 /**< right bound of interval */
   )
{
   int energy;
   int ect_j;
   int min;
   int j;

   assert(scip != NULL);
   assert(varhashmap != NULL);
   assert(var != NULL);
   assert(est < lct);

   j = (int)(size_t)SCIPhashmapGetImage(varhashmap, var);

   /* compute variable's bounds */
   ect_j = convertBoundToInt(scip, SCIPvarGetLbLocal(var)) + durations[j];

   min = MIN( lct-est, ect_j-est );

   energy =  MAX(0, MIN(min, durations[j])) * demands[j];

   return energy;
}

/** computes the energy in the interval [est,lct] of the given variable */
static
int getVarEnergy(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   SCIP_HASHMAP*         varhashmap,         /**< hashmap from variables to their indices */
   SCIP_VAR*             var,                /**< variable whose energy is reported */
   int                   est,                /**< left bound of interval */
   int                   lct                 /**< right bound of interval */
   )
{
   int energy;
   int lst_j;
   int ect_j;
   int min;
   int j;

   assert(scip != NULL);
   assert(varhashmap != NULL);
   assert(var != NULL);
   assert(est < lct);

   SCIPdebugMessage("perform energetic reasoning\n");

   j = (int)(size_t)SCIPhashmapGetImage(varhashmap, var);

   /* compute variable's bounds */
   ect_j = convertBoundToInt(scip, SCIPvarGetLbLocal(var)) + durations[j];
   lst_j = convertBoundToInt(scip, SCIPvarGetUbLocal(var));

   min = MIN3( lct-est, ect_j-est, lct-lst_j );

   energy =  MAX(0, MIN( min, durations[j] )) * demands[j];

   return energy;
}

/** computes the energy in the interval [est,lct] of all variable/jobs */
static
int computeEnergy(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   est,                /**< left bound of interval */
   int                   lct                 /**< right bound of interval */
   )
{

   int energy;
   int j;

   assert(scip != NULL);
   assert(est < lct);

   SCIPdebugMessage("perform energetic reasoning\n");

   energy = 0;

   /* loop over all jobs to compute their energetic demand in [est, lct] */
   for( j = 0; j < nvars; ++j )
   {
      SCIP_VAR* var;
      int lst_j;
      int ect_j;
      int min;

      var = vars[j];
      ect_j = convertBoundToInt(scip, SCIPvarGetLbLocal(var)) + durations[j];
      lst_j = convertBoundToInt(scip, SCIPvarGetUbLocal(var));

      min = MIN3( lct-est, ect_j-est, lct-lst_j );

      energy += MAX(0, MIN(min, durations[j])) * demands[j];
   }

   return energy;
}

/** detects whether new edges should be added to the relaxation */
static
SCIP_RETCODE performEnergeticReasoning(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_CONS*            cons,               /**< constraint which is propagated (needed to SCIPinferVar**Cons()) */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is infeasible */
   )
{
   SCIP_HASHMAP* varhashmap;

   int* ests;
   int* lcts;

   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   int est;
   int lct;

   int ntimepoints;
   int ntimepointsest;
   int ntimepointslct;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   SCIPdebugMessage("perform energetic reasoning\n");

   infeasible = FALSE;
   ntimepoints = 2*nvars;

   SCIP_CALL( SCIPallocBufferArray(scip, &ests, ntimepoints) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lcts, ntimepoints) );

   /* insert all jobs into the hashmap to find the corresponding index of the variable */
   SCIP_CALL( SCIPhashmapCreate(&varhashmap, SCIPblkmem(scip), SCIPcalcHashtableSize(nvars)) );

   /* initializing est and lct */
   for( j = 0; j < nvars; ++j )
   {
      SCIP_VAR* var;

      var = vars[j];

      assert(!SCIPhashmapExists(varhashmap, var));
      SCIP_CALL( SCIPhashmapInsert(varhashmap, var, (void*)(size_t)j) );

      lcts[2*j] = convertBoundToInt(scip, SCIPvarGetLbLocal(var)) + durations[j];
      lcts[2*j+1] = convertBoundToInt(scip, SCIPvarGetUbLocal(var)) + durations[j];

      ests[2*j] = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
      ests[2*j+1] = convertBoundToInt(scip, SCIPvarGetUbLocal(var));
   }

   /* sort the latest completion times */
   SCIPsortInt(lcts, ntimepoints);
   SCIPsortInt(ests, ntimepoints);

   j = 0;
   for( i = 0; i < ntimepoints; ++i )
   {
      if( j == 0 || lcts[i] > lcts[j-1]  )
      {
         lcts[j] = lcts[i];
         j++;
      }
   }
   ntimepointslct = j;

   j = 0;
   for( i = 0; i < ntimepoints; ++i )
   {
      if( j == 0 || ests[i] > ests[j-1]  )
      {
         ests[j] = ests[i];
         j++;
      }
   }
   ntimepointsest = j;

   for( i = 0; i < ntimepointsest && !infeasible; ++i )
   {
      est = ests[i];

      for( j = ntimepointslct-1; j >= 0 && !infeasible; --j )
      {
         int energy;
         int k;

         lct = lcts[j];

         /* only check non-empty intervals */
         if( lct <= est )
            break;

         energy = computeEnergy(scip, nvars, vars, durations, demands, est, lct);

         /* check all jobs */
         for( k = 0; k < nvars && !infeasible; k++ )
         {
            SCIP_VAR* var_k;
            int pos_k;
            int lst_k;
            int lct_k;
            int rightenergy_k;
            int energy_k;     /* energy that var_k contributes to [est,lct] */
            INFERINFO inferinfo;
            int demand_k;

            var_k = vars[k];
            pos_k = (int)(size_t)SCIPhashmapGetImage(varhashmap, var_k);

            lst_k = convertBoundToInt(scip, SCIPvarGetUbLocal(var_k));
            lct_k = lst_k + durations[pos_k];

            /* don't do anything if job does not intersect before [est;lct] */
            if( lst_k >= lct || lct_k <=est )
               continue;

            demand_k = demands[pos_k];
            energy_k = getVarEnergy(scip, durations, demands, varhashmap, var_k, est, lct) ;
            rightenergy_k = getVarRightEnergy(scip, durations, demands, varhashmap, var_k, est, lct) ;

            if( energy - energy_k > (capacity - demand_k) * (lct - est)
               && energy - energy_k + rightenergy_k > capacity * (lct - est) )
            {
               /* update lct_k */
               SCIP_Real diff;

               diff = energy - energy_k - (SCIP_Real)(capacity - demand_k) * (lct - est);
               diff /= (SCIP_Real)demand_k;

               lst_k = lct - (int)SCIPfeasCeil(scip, diff) - durations[pos_k];

               /* check if problem is already infeasible (energy overload)*/
               if( lst_k + durations[pos_k] < est )
               {
                  SCIPdebugMessage("energetic reasoning detected overload in [%d,%d]\n", est, lct);
                  inferinfo = getInferInfo(PROPRULE_4_ENERGETICREASONING, est, lct);
                  SCIP_CALL( initializeConflictAnalysisEnergeticReasoning(scip, nvars, vars, durations, NULL, inferinfo) );
                  infeasible = TRUE;
                  *initialized = TRUE;
                  *cutoff = TRUE;
               }
               else  /* problem seems feasible --> update bound */
               {
                  inferinfo = getInferInfo(PROPRULE_4_ENERGETICREASONING, est, lct);

                  SCIPdebugMessage("energetic reasoning updates var <%s>[dur=%d, dem=%d] ub from %g to %d in interval [%d,%d]\n",
                     SCIPvarGetName(var_k), durations[pos_k], demand_k,
                     SCIPvarGetUbLocal(var_k), lst_k, est, lct);
                  SCIP_CALL( SCIPinferVarUbCons(scip, var_k, (SCIP_Real)lst_k, cons, inferInfoToInt(inferinfo), TRUE, &infeasible, &tightened) );

                  if( tightened )
                     ++(*nchgbds);

                  if( infeasible )
                  {
                     SCIPdebugMessage("energetic reasoning detected infeasibility: ub-update\n");
                     SCIP_CALL( initializeConflictAnalysisEnergeticReasoning(scip, nvars, vars, durations, var_k, inferinfo) );
                     *initialized = TRUE;
                     *cutoff = TRUE;
                  }
               }
               /* recompute energy */
               energy = energy - energy_k + getVarEnergy(scip, durations, demands, varhashmap, var_k, est, lct);

            }
         }

         /* check all jobs for lb update */
         for( k = 0; k < nvars && !infeasible; ++k)
         {
            SCIP_VAR* var_k;
            int pos_k;
            int est_k;
            int ect_k;
            int energy_k;     /* energy that var_k contributes to [est,lct] */
            int leftenergy_k;
            int demand_k;
            INFERINFO inferinfo;

            var_k = vars[k];
            pos_k = (int)(size_t)SCIPhashmapGetImage(varhashmap, var_k);

            est_k = convertBoundToInt(scip, SCIPvarGetLbLocal(var_k));
            ect_k = est_k + durations[pos_k];

            /* don't do anything if job does not intersect [est;lct] */
            if( ect_k <= est || est_k >= lct )
               continue;

            demand_k = demands[pos_k];
            energy_k = getVarEnergy(scip, durations, demands, varhashmap, var_k, est, lct) ;
            leftenergy_k = getVarLeftEnergy(scip, durations, demands, varhashmap, var_k, est, lct) ;

            if( energy - energy_k > (capacity - demand_k) * (lct - est)
               && energy - energy_k + leftenergy_k > capacity * (lct - est) )
            {
               /* update est_k */
               SCIP_Real diff;

               diff = energy - energy_k - (SCIP_Real)(capacity - demand_k) * (lct - est);
               diff /= (SCIP_Real)demand_k;

               est_k = est + (int)SCIPfeasCeil(scip, diff);

               if( est_k > lct )
               {
                  SCIPdebugMessage("energetic reasoning detected overload in [%d,%d]\n", est, lct);
                  inferinfo = getInferInfo(PROPRULE_4_ENERGETICREASONING, est, lct);
                  SCIP_CALL( initializeConflictAnalysisEnergeticReasoning(scip, nvars, vars, durations, NULL, inferinfo) );
                  infeasible = TRUE;
                  *initialized = TRUE;
                  *cutoff = TRUE;
               }
               else  /* problem seems feasible --> update bound */
               {
                  inferinfo = getInferInfo(PROPRULE_4_ENERGETICREASONING, est, lct);
                  SCIPdebugMessage("energetic reasoning updates var <%s>[dur=%d, dem=%d] lb from %g to %d in interval [%d,%d]\n",
                     SCIPvarGetName(var_k),  durations[pos_k], demand_k,
                     SCIPvarGetLbLocal(var_k), est_k, est, lct);
                  SCIP_CALL( SCIPinferVarLbCons(scip, var_k, (SCIP_Real)est_k, cons, inferInfoToInt(inferinfo), TRUE, &infeasible, &tightened) );

                  if( tightened )
                     ++(*nchgbds);

                  if( infeasible )
                  {
                     SCIPdebugMessage("energetic reasoning detected infeasibility in Node %lld: lb-update\n",
                        SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

                     SCIP_CALL( initializeConflictAnalysisEnergeticReasoning(scip, nvars, vars, durations, var_k, inferinfo) );
                     *initialized = TRUE;
                     *cutoff = TRUE;
                  }
               }

               /* recompute energy */
               energy = energy - energy_k + getVarEnergy(scip, durations, demands, varhashmap, var_k, est, lct);
            }
         }

         /* go to next change in lct */
         while( j > 0 && lcts[j-1] == lct )
            --j;
      } /*lint --e{850}*/

      /* go to next change in est */
      while( i < ntimepointsest-1 && ests[i+1] == est )
         ++i;
   } /*lint --e{850}*/

   /* free hashmap */
   SCIPhashmapFree(&varhashmap);

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &lcts);
   SCIPfreeBufferArray(scip, &ests);

   return SCIP_OKAY;
}

/** repropagation of Edge finding algorithm simplified version from Petr Vilim
 *  only a small subset is reported such that energy in total and for bound change is enough
 */
static
SCIP_RETCODE analyzeShortConflictEdgeFinding(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_VAR*             infervar,           /**< variable whose bound change is to be explained */
   INFERINFO             inferinfo,          /**< inference info containing position of correct bdchgids */
   int                   inferdemand,        /**< demand of inference variable */
   int                   inferduration,      /**< duration of inference variable */
   int                   inferdiff,          /**< difference that has to be reached by edge-finder for bound update */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_Bool*            success             /**< pointer to store if we could explain the bound change */
   )
{
   SCIP_Real neededenergy;

   int* varids;          /* array of job indices that intersect the omega-interval */
   int* energies;        /* and their corresponding energies */
   int sizeenergies;

   int est_omega;
   int lct_omega;
   int delta_omega;

   int energy;
   int inferenergy;
   int j;

   assert(scip != NULL);
   assert(nvars > 0);
   assert(inferInfoGetProprule(inferinfo) == PROPRULE_3_EDGEFINDING);

   SCIPdebugMessage("repropagate edge-finding with short reasons for variable <%s>\n", SCIPvarGetName(infervar));

   *success = FALSE;
   sizeenergies = 0;

   est_omega = inferInfoGetEst(inferinfo);
   lct_omega = inferInfoGetLct(inferinfo);

   SCIP_CALL( SCIPallocBufferArray(scip, &energies, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varids, nvars) );

   /* collect the energies of all variables in [est_omega, lct_omega] */
   for( j = 0; j < nvars; ++j )
   {
      SCIP_VAR* var;
      int lb;
      int ub;

      var = vars[j];

      if( var == infervar )
         continue;

      lb = convertBoundToInt(scip, SCIPvarGetLbAtIndex(var, bdchgidx, FALSE));
      ub = convertBoundToInt(scip, SCIPvarGetUbAtIndex(var, bdchgidx, FALSE));

      /* report all jobs running in [est_omega; lct_omega] (these might be too much, but ok) */
      if( lb >= est_omega &&  ub + durations[j] <= lct_omega )
      {
         energies[sizeenergies] = durations[j] * demands[j];
         varids[sizeenergies] = j;
         sizeenergies++;
      }
   }

   SCIPsortDownIntInt(energies, varids, sizeenergies);

   delta_omega = lct_omega - est_omega;
   neededenergy = ((SCIP_Real)(capacity - inferdemand) * delta_omega) / (SCIP_Real)inferdemand;
   inferenergy = inferdemand * inferduration;

   /* report conflicting jobs until enough energy is reported */
   energy = 0;
   for( j = 0; j < sizeenergies; ++j )
   {
      SCIP_Real remaining;
      energy += energies[j];

      SCIP_CALL( SCIPaddConflictUb(scip, vars[varids[j]], bdchgidx) );
      SCIP_CALL( SCIPaddConflictLb(scip, vars[varids[j]], bdchgidx) );

      remaining = SCIPfeasCeil(scip, (SCIP_Real)energy - neededenergy);

      if( remaining >= inferdiff && energy + inferenergy > capacity * delta_omega )
      {
         *success = TRUE;
         break;
      }

#ifdef SCIP_DEBUG
      if( remaining >= inferdiff )
      {
         SCIPdebugMessage("enough energy for C-c_i\n");
      }
      if( energy + inferenergy > capacity * delta_omega )
      {
         SCIPdebugMessage("enough energy for C\n");
      }
#endif
   }

   /* free buffer array */
   SCIPfreeBufferArray(scip, &varids);
   SCIPfreeBufferArray(scip, &energies);

   if( !(*success) )
   {
      SCIPinfoMessage(scip, NULL, "could not resolve conflict from edgefinding\n");
      SCIPABORT();
   }


   return SCIP_OKAY;
}

/** repropagation of Edge finding algorithm simplified version from Petr Vilim*/
static
SCIP_RETCODE analyzeConflictEdgeFinding(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   SCIP_VAR*             infervar,           /**< variable whose bound change is to be explained */
   INFERINFO             inferinfo,          /**< inference info containing position of correct bdchgids */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_Bool*            success             /**< pointer to store if we could explain the bound change */
   )
{
   int est_omega;
   int lct_omega;
   int j;

   assert(scip != NULL);
   assert(nvars > 0);
   assert(inferInfoGetProprule(inferinfo) == PROPRULE_3_EDGEFINDING);

   SCIPdebugMessage("repropagate edge-finding variable <%s>\n", SCIPvarGetName(infervar));

   *success = FALSE;

   est_omega = inferInfoGetEst(inferinfo);
   lct_omega = inferInfoGetLct(inferinfo);

   /* collect the current lower bound of the start variables */
   for( j = 0; j < nvars; ++j )
   {
      SCIP_VAR* var;
      var = vars[j];

      if( var == infervar )
         continue;

      /* report all jobs running in [est_omega; lct_omega] (these might be too much, but ok) */
      if( convertBoundToInt(scip, SCIPvarGetLbAtIndex(var, bdchgidx, FALSE)) >= est_omega
         && convertBoundToInt(scip, SCIPvarGetUbAtIndex(var, bdchgidx, FALSE)) + durations[j] <= lct_omega )
      {
         SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx) );
         *success = TRUE;
      }
   }

   if( !(*success) )
   {
      SCIPinfoMessage(scip, NULL, "could not resolve conflict from edgefinding\n");
      SCIPABORT();
   }

   return SCIP_OKAY;
}


/** initialize conflict analysis and analyze conflict */
static
SCIP_RETCODE initializeConflictAnalysisEdgeFinding(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   SCIP_VAR*             infervar,           /**< inference variable */
   INFERINFO             inferinfo           /**< inference info */
   )
{
   SCIP_Bool success;

   SCIPdebugMessage("initialize conflict analysis\n");

   assert(inferInfoGetProprule(inferinfo) == PROPRULE_3_EDGEFINDING);

   /* conflict analysis can only be applied in solving stage and if it is turned on */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPinitConflictAnalysis(scip) );

   /* add lower and upper bound of variable which leads to the infeasibility */
   SCIP_CALL( SCIPaddConflictLb(scip, infervar, NULL) );
   SCIP_CALL( SCIPaddConflictUb(scip, infervar, NULL) );

   SCIPdebugMessage("add lower and upper bounds of variable <%s>\n", SCIPvarGetName(infervar));

   SCIP_CALL( analyzeConflictEdgeFinding(scip, nvars, vars, durations,
         infervar, inferinfo, NULL, &success) );
   assert(success);

   return SCIP_OKAY;
}

/** resolve propagation w.r.t. the cumulative condition */
static
SCIP_RETCODE respropCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_VAR*             infervar,           /**< the conflict variable whose bound change has to be resolved */
   INFERINFO             inferinfo,          /**< the user information */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   )
{
   SCIP_VAR* var;

   SCIP_Bool success;

   int inferdemand;
   int inferduration;  /* needed for upperbound resolve process */
   int j;

   assert(inferInfoGetProprule(inferinfo) != PROPRULE_INVALID);

   if( SCIPvarGetType(infervar) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(infervar) == SCIP_VARTYPE_IMPLINT )
   {
      /* get duration and demand of inference variable */
      /**@todo hashmap for variables and durations would speed this up */
      inferdemand = 0;
      inferduration = 0;

      for( j = 0; j < nvars; ++j )
      {
         var = vars[j];
         assert(var != NULL);

         if( var == infervar )
         {
            inferdemand = demands[j];
            inferduration = durations[j];
            break;
         }
      }

      SCIPdebugMessage("variable <%s> has duration = %d and demand = %d\n",
         SCIPvarGetName(infervar), inferduration, inferdemand);

      /* repropagation for core-times */
      if(  inferInfoGetProprule(inferinfo) == PROPRULE_1_CORETIMES )
      {
         int leftbound;
         int rightbound;

         if( boundtype == SCIP_BOUNDTYPE_UPPER )
         {
            SCIPdebugMessage("variable <%s> bound changed from %g to %g\n",
               SCIPvarGetName(infervar), SCIPvarGetUbAtIndex(infervar, bdchgidx, FALSE),
               SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE));

            rightbound = convertBoundToInt(scip, SCIPvarGetUbAtIndex(infervar, bdchgidx, FALSE)) + inferduration;
            leftbound = convertBoundToInt(scip, SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE)) + inferduration;

            /* old upper bound of variable itself is responsible */
            SCIP_CALL( SCIPaddConflictUb( scip, infervar, bdchgidx ) );
         }
         else
         {
            assert(boundtype == SCIP_BOUNDTYPE_LOWER);
            leftbound = convertBoundToInt(scip, SCIPvarGetLbAtIndex(infervar, bdchgidx, FALSE));
            rightbound = convertBoundToInt(scip, SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE));

            /* old lower bound of variable itself is responsible */
            SCIP_CALL( SCIPaddConflictLb( scip, infervar, bdchgidx ) );

         }
         assert(leftbound < rightbound);

         SCIP_CALL( analyzeConflictCoreTimesCumulative(scip, nvars, vars, durations, demands, capacity,
               infervar, leftbound, rightbound, inferduration, inferdemand,
               boundtype, bdchgidx, &success) );
         assert(success);
      }
      else
      {
         /* repropagation for edge-finding or energetic reasoning */
         int oldbound;
         int newbound;

         SCIPdebugMessage("repropagate edge-finder or energetic reasoning!\n");

         if( boundtype == SCIP_BOUNDTYPE_LOWER )
         {
            SCIPdebugMessage("variable <%s> lower bound changed from %g to %g\n",
               SCIPvarGetName(infervar), SCIPvarGetLbAtIndex(infervar, bdchgidx, FALSE),
               SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE));

            oldbound = convertBoundToInt(scip, SCIPvarGetLbAtIndex(infervar, bdchgidx, FALSE));
            newbound = convertBoundToInt(scip, SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE));
            assert(oldbound < newbound);

            SCIP_CALL( SCIPaddConflictLb(scip, infervar, bdchgidx) );

            /* analyze the conflict */
            if( inferInfoGetProprule(inferinfo) == PROPRULE_3_EDGEFINDING )
            {
               /* can search for small clauses if earliest start is in the interval */
               if( oldbound >= inferInfoGetEst(inferinfo) )
               {
                  int inferdiff;
                  inferdiff = newbound - inferInfoGetEst(inferinfo);
                  assert(inferdiff > 0);
                  SCIP_CALL( analyzeShortConflictEdgeFinding(scip, nvars, vars, durations, demands, capacity,
                        infervar, inferinfo, inferdemand, inferduration, inferdiff,
                        bdchgidx, &success) );
               }
               else
               {
                  SCIP_CALL( analyzeConflictEdgeFinding(scip, nvars, vars, durations,
                        infervar, inferinfo, bdchgidx, &success) );
               }
            }
            else
            {
               assert(inferInfoGetProprule(inferinfo) == PROPRULE_4_ENERGETICREASONING);

               SCIP_CALL( analyzeConflictEnergeticReasoning(scip, nvars, vars, durations,
                     infervar, inferinfo, bdchgidx, &success) );
            }
         }
         else /* now consider upper bound changes */
         {
            SCIPdebugMessage("variable <%s> upper bound changed from %g to %g\n",
               SCIPvarGetName(infervar), SCIPvarGetUbAtIndex(infervar, bdchgidx, FALSE),
               SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE));

            oldbound = convertBoundToInt(scip, SCIPvarGetUbAtIndex(infervar, bdchgidx, FALSE));
            newbound = convertBoundToInt(scip, SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE));
            assert(oldbound > newbound);

            SCIP_CALL( SCIPaddConflictUb(scip, infervar, bdchgidx) );

            /* analyze the conflict */
            if( inferInfoGetProprule(inferinfo) == PROPRULE_3_EDGEFINDING )
            {
               /* can search for small clauses if latest completion time is in the interval */
               if( oldbound + inferduration<= inferInfoGetLct(inferinfo) )
               {
                  int inferdiff;
                  inferdiff = inferInfoGetLct(inferinfo) - newbound - inferduration;
                  assert(inferdiff > 0);
                  SCIP_CALL( analyzeShortConflictEdgeFinding(scip, nvars, vars, durations, demands, capacity,
                        infervar, inferinfo, inferdemand, inferduration, inferdiff,
                        bdchgidx, &success) );
               }
               else
               {
                  SCIP_CALL( analyzeConflictEdgeFinding(scip, nvars, vars, durations,
                        infervar, inferinfo, bdchgidx, &success) );
               }
            }
            else /* upper bound conflict analysis for energetic reasoning */
            {
               assert(inferInfoGetProprule(inferinfo) == PROPRULE_4_ENERGETICREASONING);

               SCIP_CALL( analyzeConflictEnergeticReasoning(scip, nvars, vars, durations,
                     infervar, inferinfo, bdchgidx, &success) );
            }
         }
         assert(success);
      }
   }
   else
   {
      /* repropagation for binary variables set to zero; inferinfo == position in array and excluded timepoint */
      SCIP_VAR* intvar;
      SCIP_VAR** binvars;

      int pos;
      int nbinvars;

      assert(SCIPvarGetType(infervar) == SCIP_VARTYPE_BINARY);
      assert(inferInfoGetProprule(inferinfo) == PROPRULE_2_CORETIMEHOLES);

      intvar = NULL;
      inferdemand = 0;

      pos = inferInfoGetEst(inferinfo);
      assert(pos >= 0);

      /* get demand and integer variable of given inference variable */
      for( j = 0; j < nvars; ++j )
      {
         var = vars[j];
         assert(var != NULL);

         SCIP_CALL( SCIPgetBinvarsLinking(scip, SCIPgetConsLinking(scip, var), &binvars, &nbinvars) );

         if( binvars[pos] == infervar )
         {
            intvar = var;
            inferdemand = demands[j];
            break;
         }
      }
      assert(intvar != NULL);
      assert(inferdemand > 0);

      SCIP_CALL( analyzeConflictCoreTimesBinvarsCumulative(scip, nvars, vars, durations, demands, capacity,
            infervar, intvar, inferInfoGetLct(inferinfo), inferdemand, bdchgidx, &success) );
   }

   if( success )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/** computes a new earliest starting time of the job in 'respleaf' due to the energy consumption and stores the
 *  responsible interval bounds in *est_omega and *lct_omega
 */
static
int computeNewLstOmegaset(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_HASHMAP*         varhashmap,         /**< hashmap variable -> index */
   TLTREENODE*           respleaf,           /**< theta tree leaf whose variables lower bound will be updated */
   TLTREENODE**          omegaset,           /**< set of lambda nodes who determine new lower bound */
   int                   nelements,          /**< number of elements in omegaset */
   int                   lct_j,              /**< latest completion time from all nodes in theta lambda tree */
   int                   makespan,           /**< upper bound on lct of all vars */
   int*                  est_omega,          /**< pointer to store est of set omega */
   int*                  lct_omega           /**< pointer to store lct of set omega */
   )
{
   int newest;
   int newlst;
   int tmp;
   int j;
   int pos;
   int energy;

   int demand;
   int duration;

   assert(scip != NULL);

   energy = 0;
   newest = 0;
   *est_omega = INT_MAX;
   *lct_omega = 0;

   /* get position of responsible variable */
   pos = (int)(size_t)SCIPhashmapGetImage(varhashmap, respleaf->var);
   assert(pos < nvars);

   demand = demands[pos];
   duration = durations[pos];
   assert(demand > 0);

   for( j = 0; j < nelements; ++j )
   {
      int idx;

      assert(omegaset[j]->inTheta);
      assert(omegaset[j]->var != NULL);
      idx = (int)(size_t)SCIPhashmapGetImage(varhashmap, omegaset[j]->var);
      assert(idx < nvars);

      tmp = (int) (makespan - SCIPvarGetUbLocal(omegaset[j]->var) - durations[idx]);
      *est_omega = MIN(*est_omega, tmp);
      tmp = (int) (makespan - SCIPvarGetLbLocal(omegaset[j]->var));
      *lct_omega = MAX(*lct_omega, tmp);

      assert(durations[idx] * demands[idx] == omegaset[j]->energy);
      assert(*lct_omega <= lct_j);
      energy += omegaset[j]->energy;
   }

   /* update est if enough energy */
   if( energy > (capacity - demand) * (*lct_omega - *est_omega) )
   {
      if( energy + demand * duration > capacity * (*lct_omega - *est_omega) )
      {
         newest = (int)SCIPfeasCeil(scip, (energy - (SCIP_Real)(capacity - demand) * (*lct_omega - *est_omega)) / (SCIP_Real)demand);
         newest += *est_omega;
      }

      assert(energy + demand * duration > capacity * (*lct_omega - *est_omega));

      /* recompute original values using 'makespan' */
      tmp = makespan - *est_omega;
      *est_omega =  makespan - *lct_omega;
      *lct_omega = tmp;
   }

   newlst = makespan - newest - durations[pos];

   return newlst;
}

/** computes a new latest starting time of the job in 'respleaf' due to the energy consumption and stores the
 *  responsible interval bounds in *est_omega and *lct_omega
 */
static
int computeNewEstOmegaset(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_HASHMAP*         varhashmap,         /**< hashmap variable -> index */
   TLTREENODE*           respleaf,           /**< theta tree leaf whose variables lower bound will be updated */
   TLTREENODE**          omegaset,           /**< set of lambda nodes who determine new lower bound */
   int                   nelements,          /**< number of elements in omegaset */
   int                   lct_j,              /**< latest completion time from all nodes in theta lambda tree */
   int*                  est_omega,          /**< pointer to store est of set omega */
   int*                  lct_omega           /**< pointer to store lct of set omega */
   )
{
   int newest;
   int j;
   int pos;
   int energy;

   int demand;
   int duration;

   assert(scip != NULL);

   energy = 0;
   newest = 0;
   *est_omega = INT_MAX;
   *lct_omega = 0;

   /* get position of responsible variable */
   pos = (int)(size_t)SCIPhashmapGetImage(varhashmap, respleaf->var);
   assert(pos < nvars);

   demand = demands[pos];
   duration = durations[pos];
   assert(demand > 0);

   for( j = 0; j < nelements; ++j )
   {
      int idx ;
      int tmp;

      assert(omegaset[j]->var != NULL);
      assert(omegaset[j]->inTheta);
      idx = (int)(size_t)SCIPhashmapGetImage(varhashmap, omegaset[j]->var);
      assert(idx < nvars);

      tmp = convertBoundToInt(scip, SCIPvarGetLbLocal(omegaset[j]->var));
      *est_omega = MIN(*est_omega, tmp);
      tmp = convertBoundToInt(scip, SCIPvarGetUbLocal(omegaset[j]->var)) + durations[idx];
      *lct_omega = MAX(*lct_omega, tmp);

      assert(durations[idx] * demands[idx] == omegaset[j]->energy);
      assert(*lct_omega <= lct_j);

      energy += omegaset[j]->energy;
   }

   if( energy >  (capacity - demand) * (*lct_omega - *est_omega) )
   {
      if( energy + demand * duration > capacity * (*lct_omega - *est_omega) )
      {
         newest =  (int)SCIPfeasCeil(scip, (energy - (SCIP_Real)(capacity - demand) * (*lct_omega - *est_omega)) / (SCIP_Real)demand);
         newest += (*est_omega);
      }
      assert(energy + demand * duration > capacity * ( (*lct_omega) - (*est_omega) ));
   }

   return newest;
}

/** detects whether new edges should be added to the relaxation */
static
SCIP_RETCODE performEdgeFindingDetection(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             forward,            /**< shall lower bounds be updated? */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_CONS*            cons,               /**< constraint which is propagated (needed to SCIPinferVar**Cons()) */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is infeasible */
   )
{
   TLTREENODE** nodes;
   TLTREE* tltree;
   SCIP_HASHMAP* varhashmap;

   int* lcts;
   int* lct_ids;

   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   int makespan;             /* needed if we want to make the lct-updates instead of est */
   int j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff !=  NULL);
   assert(*cutoff ==  FALSE);
   assert(initialized !=  NULL);
   assert(*initialized ==  FALSE);


   SCIPdebugMessage("perform edge-finding for cumulative condition of constraint <%s>\n", SCIPconsGetName(cons));

   infeasible = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &lcts, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lct_ids, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes, nvars) );

   /* insert all jobs into the hashmap to find the corresponding index of the variable */
   SCIP_CALL( SCIPhashmapCreate(&varhashmap, SCIPblkmem(scip), SCIPcalcHashtableSize(nvars)) );

   /* compute makespan */
   makespan = 0;
   if( !forward )
   {
      for( j = 0; j < nvars; ++j )
      {
         int tmp;

         tmp = convertBoundToInt(scip, SCIPvarGetUbLocal(vars[j])) + durations[j];
         makespan = MAX(makespan, tmp);
      }
   }

   /* initializing latest completion times and tree leaves */
   for( j = 0; j < nvars; ++j )
   {
      SCIP_VAR* var;
      SCIP_Real est;
      int energy;

      var = vars[j];
      assert(var != NULL);
      assert(!SCIPhashmapExists(varhashmap, var));

      SCIP_CALL( SCIPhashmapInsert(varhashmap, var, (void*)(size_t)j) );

      if( forward )
      {
         lcts[j] = convertBoundToInt(scip, SCIPvarGetUbLocal(var)) + durations[j];
         lct_ids[j] = j;

         est = convertBoundToInt(scip, SCIPvarGetLbLocal(var) ) + j / (2.0 * nvars);
         energy = demands[j] * durations[j];
      }
      else
      {
         lcts[j] = makespan - convertBoundToInt(scip, SCIPvarGetLbLocal(var));
         lct_ids[j] = j;

         est = makespan - convertBoundToInt(scip, SCIPvarGetUbLocal(var)) - durations[j] + j / (2.0 * nvars);
         energy = demands[j] * durations[j];
      }

      SCIP_CALL( tltreeCreateThetaLeaf(scip, &(nodes[j]), var, est, energy, capacity * (int)(est + 0.01) + energy) );
   }

   /* sort the latest completion times */
   SCIPsortIntInt(lcts, lct_ids, nvars);

   SCIP_CALL( tltreeCreateTree(scip, &tltree, nodes, lct_ids, nvars) );

   /* iterate over all jobs in non-decreasing order of lct_j */
   for( j = nvars-1; !infeasible && j >= 0; --j )
   {
      while( !infeasible && tltreeGetEnvelopTL(tltree) > capacity * lcts[j] )
      {
         TLTREENODE** omegaset;
         TLTREENODE* respleaf;

         int nelements;
         INFERINFO inferinfo;

         int est_omega;
         int lct_omega;

         int pos;
         int duration_pos;

         /* find out which variable var_i from Lambda is responsible */
         respleaf = tltreeFindResponsibleLeaf(tltree);

         assert(respleaf != NULL);
         assert(respleaf->left == NULL);
         assert(respleaf->right == NULL);
         assert(respleaf->var != NULL);
         assert(respleaf->energyL > 0);

         /* get position of responsible variable */
         pos = (int)(size_t)SCIPhashmapGetImage(varhashmap, respleaf->var);
         assert(pos < nvars);

         duration_pos = durations[pos];

         if( respleaf->value + duration_pos/*respleaf->energyL*/ >= lcts[j] )
         {
            SCIP_CALL( tltreeDeleteLeaf(scip, tltree, respleaf) );
            continue;
         }

         /* compute omega set */
         SCIP_CALL( SCIPallocBufferArray(scip, &omegaset, nvars) );

         /* This was the previous line. In my case I got invalid write which means the array was to small
          * SCIP_CALL( SCIPallocBufferArray(scip, &omegaset, nvars - j) );
          * ????????????????????????????????
          */

         /* get omega set from tltree */
         SCIP_CALL( tltreeReportOmegaSet(scip, tltree, &omegaset, &nelements) );

         assert(nelements > 0);
         assert(nelements < nvars);

         /* compute new earliest starting time */
         if( forward )
         {
            int newest;
            newest = computeNewEstOmegaset(scip, nvars, durations, demands, capacity,
               varhashmap, respleaf, omegaset, nelements, lcts[j], &est_omega, &lct_omega);

            inferinfo = getInferInfo(PROPRULE_3_EDGEFINDING, est_omega, lct_omega);

            /* update variable's lower bound */
            SCIP_CALL( SCIPinferVarLbCons(scip, respleaf->var, (SCIP_Real)newest,
                  cons, inferInfoToInt(inferinfo), TRUE, &infeasible, &tightened) );
         }
         else
         {
            int newlst;
            newlst = computeNewLstOmegaset(scip,  nvars, durations, demands, capacity,
               varhashmap, respleaf, omegaset, nelements, lcts[j], makespan, &est_omega, &lct_omega);

            inferinfo = getInferInfo(PROPRULE_3_EDGEFINDING, est_omega, lct_omega);

            /* update variable's upper bound */
            SCIP_CALL( SCIPinferVarUbCons(scip, respleaf->var, (SCIP_Real)newlst,
                  cons, inferInfoToInt(inferinfo), TRUE, &infeasible, &tightened) );
         }

         /* if node can be cutoff, start conflict analysis */
         if( infeasible )
         {
            SCIP_CALL( initializeConflictAnalysisEdgeFinding(scip, nvars, vars, durations,
                  respleaf->var, inferinfo) );
            *initialized = TRUE;
            *cutoff = TRUE;
         }

         /* count number of tightened bounds */
         if( tightened )
            ++(*nchgbds);

         /* free omegaset array */
         SCIPfreeBufferArray(scip, &omegaset);

         /* delete responsible leaf from lambda */
         SCIP_CALL( tltreeDeleteLeaf(scip, tltree, respleaf) );
      }

      /* change set of job j from Theta to Lambda */
      SCIP_CALL( tltreeTransformLeafTtoL(scip, tltree, nodes[lct_ids[j]]) );
   }

   /* free theta tree */
   SCIP_CALL( freeTltree(scip, &tltree) );

   for( j = 0; j < nvars; ++j )
   {
      if( nodes[j] != NULL )
      {
         SCIP_CALL( freeTltreeLeaf(scip, &(nodes[j])) );
      }
   }

   /* free hashmap */
   SCIPhashmapFree(&varhashmap);

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &nodes);
   SCIPfreeBufferArray(scip, &lct_ids);
   SCIPfreeBufferArray(scip, &lcts);

   return SCIP_OKAY;
}

/** checks whether the instance is infeasible due to overload,
 *  see Vilim: CPAIOR 2009: Max Energy Filtering Algorithm for Discrete Cumulative Resources
 */
static
SCIP_RETCODE checkOverload(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_CONS*            cons,               /**< constraint which is propagated */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is infeasible */
   )
{
   THETATREENODE** nodes;
   THETATREE* thetatree;

   int* lcts;
   int* lct_ids;

   int j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(initialized != NULL);
   assert(cutoff != NULL);

   SCIPdebugMessage("check/propagate overload of cumulative condition of constraint <%s>\n", SCIPconsGetName(cons));

   SCIP_CALL( SCIPallocBufferArray(scip, &lcts, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lct_ids, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes, nvars) );

   /* initializing latest completion times */
   for( j = 0; j < nvars; ++j )
   {
      SCIP_VAR* var;
      SCIP_Real est;
      int energy;

      var = vars[j];
      lcts[j] = convertBoundToInt(scip, SCIPvarGetUbLocal(var)) + durations[j];
      lct_ids[j] = j;

      est = convertBoundToInt(scip, SCIPvarGetLbLocal(var)) + j / (2.0 * nvars);
      energy = demands[j] * durations[j];

      SCIP_CALL( thetatreeCreateLeaf(scip, &(nodes[j]), var, est,
            energy, capacity * (int)(est + 0.01) + energy ) );
   }

   /* sort the latest completion times */
   SCIPsortIntInt(lcts, lct_ids, nvars);

   SCIP_CALL( createThetaTree(scip, &thetatree) );

   /* iterate over all jobs in non-decreasing order of lct_j */
   for( j = 0; j < nvars && !(*cutoff); ++j )
   {
      SCIP_Bool inserted;
      SCIP_CALL( thetatreeInsertLeaf(scip, thetatree, nodes[lct_ids[j]], &inserted) );
      assert(inserted);

      if( thetaTreeGetEnvelop(thetatree) > capacity * lcts[j] )
      {
         /*@todo: start conflict analysis, compute conflicting set */
         SCIPdebugMessage("Overload detected! Node can be cut off @todo: start conflict analysis\n");
         *initialized = FALSE;
         *cutoff = TRUE;
      }
   }

   /* free theta tree */
   SCIP_CALL( freeThetaTree(scip, &thetatree) );

   for( j = 0; j < nvars; ++j )
   {
      if( nodes[j] != NULL )
      {
         SCIP_CALL( freeThetaTreeLeaf(scip, &(nodes[j])) );
      }
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &nodes);
   SCIPfreeBufferArray(scip, &lct_ids);
   SCIPfreeBufferArray(scip, &lcts);

   return SCIP_OKAY;
}


/** fix integer variable to lower bound if the rounding locks and the object coefficient are in favor of that */
static
SCIP_RETCODE fixIntegerVariableLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< integer variable to fix */
   int*                  nchgbds             /**< pointer to store the number changed variable bounds */
   )
{
   SCIP_Real objval;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   assert(!SCIPinRepropagation(scip));

   /* if SCIP is in probing mode we cannot perform this dual reductions since this dual reduction would end in an
    * implication which can lead to cutoff the optimal solution
    */
   if( SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* rounding the variable to the lower bound is only a feasible dual reduction if the cumulative constraint
    * handler is the only one locking that variable down
    */
   assert( SCIPvarGetNLocksDown(var) >= 1 );
   if( SCIPvarGetNLocksDown(var) > 1 )
      return SCIP_OKAY;

   objval = SCIPvarGetObj(var);

   /* rounding the integer variable down is only a valid dual reduction if the object coefficient is zero or positive
    * (the transformed problem is always a minimization problem)
    */
   if( SCIPisNegative(scip, objval) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPfixVar(scip, var, SCIPvarGetLbLocal(var), &infeasible, &tightened) );
   assert(!infeasible);

   if( tightened )
      (*nchgbds)++;

   return SCIP_OKAY;
}

/** fix integer variable to upper bound if the rounding locks and the object coefficient are in favor of that */
static
SCIP_RETCODE fixIntegerVariableUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< integer variable to fix */
   int*                  nchgbds             /**< pointer to store the number changed variable bounds */
   )
{
   SCIP_Real objval;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   assert(!SCIPinRepropagation(scip));

   /* if SCIP is in probing mode we cannot perform this dual reductions since this dual reduction would end in an
    * implication which can lead to an infeasible cutoff (optimal solution)
    */
   if( SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* rounding the variable to the upper bound is only a feasible dual reduction if the cumulative constraint
    * handler is the only one locking that variable up
    */
   assert(SCIPvarGetNLocksUp(var) >= 1);
   if( SCIPvarGetNLocksUp(var) > 1 )
      return SCIP_OKAY;

   objval = SCIPvarGetObj(var);

   /* rounding the integer variable up is only a valid dual reduction if the object coefficient is zero or negative
    * (the transformed problem is always a minimization problem)
    */
   if( SCIPisPositive(scip, objval) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPfixVar(scip, var, SCIPvarGetUbLocal(var), &infeasible, &tightened) );
   assert(!infeasible);

   if( tightened )
      (*nchgbds)++;

   return SCIP_OKAY;
}

/** depending on the objective coefficient of the integer variable and the rounding locks, we might can fix the integer
 *  variable (dual reduction)
 */
static
SCIP_RETCODE fixIntegerVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< integer variable to fix */
   int*                  nchgbds             /**< pointer to store the number changed variable bounds */
   )
{
   SCIP_Real objval;
   SCIP_Real fixvalue;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   objval = SCIPvarGetObj(var);
   fixvalue = SCIP_INVALID;

   /* if SCIP is in probing mode or during repropagation we cannot perform this dual reductions since this dual
    * reduction would end in an implication which can lead to cutoff the optimal solution
    */
   if( SCIPinProbing(scip) || SCIPinRepropagation(scip) )
      return SCIP_OKAY;

   assert(SCIPvarGetNLocksDown(var) >= 1);
   assert(SCIPvarGetNLocksUp(var) >= 1);

   if( SCIPisZero(scip, objval) )
   {
      /* the integer start time variable has a zero objective value; if only the cumulative constraint handler has a
       * problem with rounding it down or up; therefore, rounding the integer down or up up is feasible dual reduction
       */
      if( SCIPvarGetNLocksDown(var) == 1 )
         fixvalue = SCIPvarGetLbLocal(var);
      else if( SCIPvarGetNLocksUp(var) == 1 )
         fixvalue = SCIPvarGetUbLocal(var);
      else
         return SCIP_OKAY;
   }
   else if( SCIPisNegative(scip, objval) && SCIPvarGetNLocksUp(var) == 1 )
   {
      /* the integer start time variable has a negative objective value and only the cumulative constraint handler has a
       * problem with rounding it up; therefore rounding it to the upper bound is the best thing we can do
       */
      fixvalue = SCIPvarGetUbLocal(var);
   }
   else if( SCIPisPositive(scip, objval) && SCIPvarGetNLocksDown(var) == 1 )
   {
      /* the integer start time variable has a positive objective value and only the cumulative constraint handler has a
       * problem with rounding it down; therefore rounding it to the lower bound is the best thing we can do
       */
      fixvalue = SCIPvarGetLbLocal(var);
   }
   else
      return SCIP_OKAY;

   /* the integer start time variable has a positive objective value and only the cumulative constraint handler has a
    * problem with rounding it down; therefore rounding it to the lower bound is the best thing we can do
    */
   assert(fixvalue < SCIP_INVALID);
   SCIP_CALL( SCIPfixVar(scip, var, fixvalue, &infeasible, &tightened) );
   assert(!infeasible);

   if( tightened )
      (*nchgbds)++;

   return SCIP_OKAY;
}

/** deletes coefficient at given position from constraint data */
static
SCIP_RETCODE deletePos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(pos < consdata->nvars);

   /* remove the rounding locks for the deleted variable */
   SCIP_CALL( unlockRounding(scip, cons, consdata->vars[pos]) );

   if( consdata->linkingconss != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &consdata->linkingconss[pos]) );
   }

   SCIPdebugMessage("remove variable <%s> from cumulative constraint <%s>\n",
      SCIPvarGetName(consdata->vars[pos]), SCIPconsGetName(cons));

   if( pos != consdata->nvars - 1 )
   {
      consdata->vars[pos] = consdata->vars[consdata->nvars-1];
      consdata->demands[pos] = consdata->demands[consdata->nvars-1];
      consdata->durations[pos] = consdata->durations[consdata->nvars-1];

      if( consdata->linkingconss != NULL )
      {
         consdata->linkingconss[pos]= consdata->linkingconss[consdata->nvars-1];
      }
   }

   consdata->nvars--;

   return SCIP_OKAY;
}

/** remove variable from the constraint and adjust permutation array */
static
SCIP_RETCODE removeVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   int                   pos,                /**< position to remove */
   int                   nvars,              /**< number of variables */
   int*                  perm                /**< permutation array to adjust */
   )
{
   /* delete variable at the given position */
   SCIP_CALL( deletePos(scip, cons, pos) );

   SCIPdebugMessage("perm[%d] = %d\n", pos, perm[nvars-1]);
   SCIPdebugMessage("perm[%d] = %d\n", perm[nvars-1], pos);

   /* fix permutation array */
   perm[pos] = perm[nvars-1];
   perm[perm[nvars-1]] = pos;

   return SCIP_OKAY;
}

/** remove irrelevant jobs; a job is irrelevant if:
 *
 *  (1) Let the jobs be non-decreasing ordered by there earlier start time and j1 and j2 first two jobs of that ordering
 *
 *     (a) if the latest completion time (lct) of job j1 is less than or equal to the earlier start time (est) of job
 *         j2, then job j1 can be scheduled at any time of the feasible time window without interfering with other jobs
 *         => j1 can be removed from the cumulative constraint
 *
 *     (b) if the earliest completion time (ect) of job j1 is less than or equal to the earliest start time (est) of job
 *         j2 and fixing the start time variable of job j1 to the lower bound is a feasible dual reduction
 *         => j1 can be removed from the cumulative constraint and fixed to its earliest start time
 *
 *  (2) Let the jobs be non-increasing ordered by there latest completion time and j1 and j2 first two jobs of that
 *      ordering
 *
 *     (a) if the earliest start time (est) of job j1 is greater than or equal to the latest completion time (lct) of
 *         job j2, then job j1 can be scheduled at any time of the feasible time window without interfering with other
 *         jobs
 *         => j1 can be removed from the cumulative constraint
 *
 *     (b) if the latest start time (lst) of job j1 is greater than or equal to the latest completion time (lct) of job
 *         j2 and fixing the start time variable of job j1 to the upper bound is a feasible dual reduction
 *         => j1 can be removed from the cumulative constraint and fixed to its latest start time
 */
static
SCIP_RETCODE removeIrrelevantJobs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   int*                  nchgbds             /**< pointer to store the number changed variable bounds */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
#if 0
   int* demands;
#endif
   int* times;
   int* indices;
   int* perm;

#if 0
   int capacity;
#endif
   int nvars;
   int idx;
   int v;

   SCIPdebugMessage("remove irrelevant jobs of cumulative constraint <%s>\n", SCIPconsGetName(cons));

   SCIPdebug( SCIP_CALL(SCIPprintCons(scip, cons, NULL) ) );

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
#if 0
   demands = consdata->demands;
   capacity = consdata->capacity;
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &times, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &indices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nvars) );

   for( v = 0; v < nvars; ++v )
   {
      times[v] = convertBoundToInt(scip, SCIPvarGetLbGlobal(consdata->vars[v]));
      indices[v] = v;
      perm[v] = v;
   }

   /* sort the indices w.r.t. to the earliest starting time */
   SCIPsortIntInt(times, indices, nvars);

   for( v = 1; v < consdata->nvars; ++v )
   {
      int duration;
      int est;
      int ect;
      int lct;

      idx = perm[indices[v-1]];
      var = consdata->vars[idx];
      assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);

      duration = consdata->durations[idx];
      ect = times[v-1] + duration;
      lct = convertBoundToInt(scip, SCIPvarGetUbGlobal(var)) + duration;
      est = times[v];
      assert(ect <= lct);
      assert(times[v-1] == convertBoundToInt(scip, SCIPvarGetLbGlobal(var)));

      if( lct <= est )
      {
         /* case (1a): job can be remove from the constraint and a dual reduction on the integer start time variable can
          * be tried
          */

         SCIPdebugMessage("variable <%s> is irrelevant\n", SCIPvarGetName(consdata->vars[idx]));

         /* fix integer start time variable if possible */
         if( SCIPconsIsChecked(cons) )
         {
            SCIP_CALL( fixIntegerVariable(scip, var, nchgbds) );
         }

         /* remove the variable and adjust the permutation array */
         SCIP_CALL( removeVariable(scip, cons, idx, consdata->nvars, perm) );
      }
      else if( ect <= est )
      {
         /* case (1b): job can be remove form the constraint only if the integer start time variable can be fixed to its
          * lower bound; fixing the integer start time variable to it lower means that afterwards case (1a) can be
          * applied
          */

         /* fix integer start time variable if possible to it lower bound */
         if( SCIPconsIsChecked(cons) )
         {
            SCIP_CALL( fixIntegerVariableLb(scip, var, nchgbds) );
         }

         if( SCIPvarGetLbGlobal(var) + 0.5 > SCIPvarGetUbGlobal(var) )
         {
            SCIPdebugMessage("variable <%s> is irrelevant\n", SCIPvarGetName(consdata->vars[idx]));

            SCIP_CALL( removeVariable(scip, cons, idx, consdata->nvars, perm) );
         }
         else
            break;
      }
      else
      {
#if 0
         int cumudemand;
         int i;
         int t;

         /* case (1c): job can be */

         /* check for each time point in the interval [est,ect) if the capacity of the resource cannot be violated */

         cumudemand = demands[idx];

         for( t = est; t < ect; ++t )
         {
            for( i = v; v < consdata->nvars && cumudemand <= capacity; ++v )
            {
               if( time[i] > t )
                  break;

               cumudemand += demands[perm[i]];
            }

            if( cumudemand > capacity )
               break;
         }

         if( t == ect )
         {
            /* fix integer start time variable if possible to it lower bound */
            if( SCIPconsIsChecked(cons) )
            {
               SCIP_CALL( fixIntegerVariableLb(scip, var, nchgbds) );
            }

            if( SCIPvarGetLbGlobal(var) + 0.5 > SCIPvarGetUbGlobal(var) )
            {
               SCIP_VAR* newvar;
               char name[SCIP_MAXSTRLEN];

               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s'", SCIPvarGetName(var));

               /* remove locks from the variable */
               SCIP_CALL( unlockRounding(scip, cons, var) );

               SCIP_CALL( SCIPcreateVar(scip, &newvar, name, est, est, 0.0, SCIP_VARTYPE_INTEGER,
                     SCIPvarIsInitial(var), SCIPvarIsRemovable(var), NULL, NULL, NULL, NULL, NULL) );
               SCIP_CALL( SCIPaddVar(scip, newvar) );
               SCIP_CALL( lockRounding(scip, cons, newvar) );

               /* replace with new variable */
               consdata->vars[idx] = newvar;
               consdata->durations[idx] = lct - est;

               SCIP_CALL( SCIPreleaseVar(scip, &newvar) );
            }
         }
#endif
         break;
      }
   }

   nvars = consdata->nvars;

   for( v = 0; v < nvars; ++v )
   {
      times[v] = convertBoundToInt(scip, SCIPvarGetUbGlobal(consdata->vars[v]) + consdata->durations[v]);
      indices[v] = v;
      perm[v] = v;
   }

   /* sort the indices w.r.t. to the latest completion time */
   SCIPsortDownIntInt(times, indices, nvars);

   for( v = 1; v < consdata->nvars; ++v )
   {
      int                                                 lct;
      int lst;
      int est;

      idx = perm[indices[v-1]];
      var = consdata->vars[idx];
      assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);

      lct = times[v];
      lst = times[v-1] - consdata->durations[idx];
      est = convertBoundToInt(scip, SCIPvarGetLbGlobal(var));
      assert(est <= lst);
      assert(lst == convertBoundToInt(scip, SCIPvarGetUbGlobal(var)));

      if( est >= lct )
      {
         /* case (2a): job can be remove form the constraint and a dual reduction on the integer start time variable can
          * be tried
          */

         SCIPdebugMessage("variable <%s> is irrelevant\n", SCIPvarGetName(consdata->vars[idx]));

         /* fix integer start time variable if possible */
         if( SCIPconsIsChecked(cons) )
         {
            SCIP_CALL( fixIntegerVariable(scip, var, nchgbds) );
         }

         SCIP_CALL( removeVariable(scip, cons, idx, consdata->nvars, perm) );
      }
      else if( lst >= lct )
      {
         /* case (2b): job can be remove form the constraint only if the integer start time variable can
          * be fixed to its upper bound
          */

         /* fix integer start time variable if possible to its upper bound */
         if( SCIPconsIsChecked(cons) )
         {
            SCIP_CALL( fixIntegerVariableUb(scip, var, nchgbds) );
         }

         if( SCIPvarGetLbGlobal(var) + 0.5 > SCIPvarGetUbGlobal(var) )
         {
            SCIPdebugMessage("variable <%s> is irrelevant\n", SCIPvarGetName(consdata->vars[idx]));

            SCIP_CALL( removeVariable(scip, cons, idx, consdata->nvars, perm) );
         }
         else
            break;
      }
      else
         break;
   }
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &indices);
   SCIPfreeBufferArray(scip, &times);

   SCIPdebugMessage("constraint <%s> contains of %d jobs\n", SCIPconsGetName(cons), consdata->nvars);

   return SCIP_OKAY;
}

/** divides demands by their greatest common divisor and divides capacity by the same value, rounding down the result */
static
void normalizeDemands(
   SCIP_CONS*            cons,               /**< knapsack constraint */
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients */
   int*                  nchgsides           /**< pointer to count number of side changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Longint gcd;
   int capacity;
   int mindemand;
   int nvars;
   int v;

   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(!SCIPconsIsModifiable(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   capacity = consdata->capacity;

   if( capacity == 1 )
      return;

   nvars = consdata->nvars;
   assert(nvars >= 1);

   /**@todo sort items w.r.t. the demands, because we can stop earlier if the smaller weights are evaluated first */

   gcd = (SCIP_Longint)consdata->demands[nvars-1];
   mindemand = 2*consdata->demands[nvars-1];

   for( v = nvars-2; v >= 0 && (gcd >= 2 || mindemand > capacity); --v )
   {
      gcd = SCIPcalcGreComDiv(gcd, (SCIP_Longint)consdata->demands[v]);
      mindemand = MIN(mindemand, 2*consdata->demands[v]);
   }

   if( mindemand > capacity )
   {
      SCIPdebugMessage("cumulative constraint <%s>: change to unary demands\n", SCIPconsGetName(cons));

      for( v = 0; v < nvars; ++v )
         consdata->demands[v] = 1;

      consdata->capacity = 1;

      (*nchgcoefs) += nvars;
      (*nchgsides)++;
   }
   else if( gcd >= 2 )
   {
      SCIPdebugMessage("cumulative constraint <%s>: dividing demands by %"SCIP_LONGINT_FORMAT"\n", SCIPconsGetName(cons), gcd);

      for( v = 0; v < nvars; ++v )
         consdata->demands[v] /= gcd;

      consdata->capacity /= gcd;

      (*nchgcoefs) += nvars;
      (*nchgsides)++;
   }
}

/** check if cumulative constraint is independently of all other constraints */
static
SCIP_Bool isConsIndependently(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< cumulative constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   vars = consdata->vars;

   /* check if the cumulative constraint has the only locks on the involved variables */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;

      var = vars[v];
      assert(var != NULL);

      if( SCIPvarGetNLocksDown(var) > 1 || SCIPvarGetNLocksUp(var) > 1 )
         return FALSE;
   }

   return TRUE;
}


/** in case the cumulative constraint is independent of every else, solve the cumulative problem and apply the fixings
 *  (dual reductions)
 */
static
SCIP_RETCODE dualPresolving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   SCIP_Longint          maxnodes,           /**< maximum nodes to solve an independent cumulative constraint (-1: unlimited) */
   int*                  nchgbds,            /**< pointer to store the number changed variable bounds */
   int*                  nfixedvars,         /**< pointer to count number of fixings */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints  */
   SCIP_Bool*            cutoff,             /**< pointer to store if the constraint is infeasible */
   SCIP_Bool*            unbounded           /**< pointer to store if the constraint is unbounded */
   )
{
   SCIP* subscip;
   SCIP_CONSDATA* consdata;
   SCIP_HASHMAP* varmapfw;
   SCIP_CONS* targetcons;
   SCIP_VAR** vars;
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   SCIP_Bool succeed;
   char probname[SCIP_MAXSTRLEN];
   int nvars;
   int v;

   assert(!SCIPconsIsModifiable(cons));
   assert(SCIPgetNConss(scip) > 0);

   /* if the cumulative constraint is the only constraint do nothing */
   if( SCIPgetNConss(scip) == 1 )
      return SCIP_OKAY;

   /* constraints for which the check flag is set to FALSE, did not contribute to the lock numbers; therefore, we cannot
    * use the locks to decide for a dual reduction using this constraint;
    */
   if( !SCIPconsIsChecked(cons) )
      return SCIP_OKAY;

   /* check if constraint is independently */
   if( !isConsIndependently(scip, cons) )
      return SCIP_OKAY;

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
   if( timelimit <= 0.0 || memorylimit <= 0.0 )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   vars = consdata->vars;

   SCIPdebugMessage("the cumulative constraint <%s> is independent from rest of the problem\n", SCIPconsGetName(cons));
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );

   /* copy all plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* get name of the original problem and add the string "_cumulative" */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_cumulative", SCIPgetProbName(scip));

   /* create the subproblem */
   SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* copy cumulative constraint */
   SCIP_CALL( SCIPgetConsCopy(scip, subscip, cons, &targetcons, SCIPconsGetHdlr(cons), varmapfw, NULL, NULL,
         FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, &succeed) );

   if( succeed )
   {
      /* add constraint to subscip */
      SCIP_CALL( SCIPaddCons(subscip, targetcons) );

      /* do not abort subproblem on CTRL-C */
      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

      /* disable output to console */
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

      /* set limits for the subproblem */
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", maxnodes) );
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

      /* forbid recursive call of heuristics and separators solving subMIPs */
      SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

      /* solve single cumulative constraint by branch and bound */
      SCIP_CALL( SCIPsolve(subscip) );

      /* evaluated solution status */
      switch( SCIPgetStatus(subscip) )
      {
      case SCIP_STATUS_INFORUNBD:
      case SCIP_STATUS_INFEASIBLE:
         *cutoff = TRUE;
         break;
      case SCIP_STATUS_UNBOUNDED:
         *unbounded = TRUE;
         break;
      case SCIP_STATUS_OPTIMAL:
      {
         /* copy optimal as dual reduction into the original SCIP instance */
         SCIP_SOL* sol;

         sol = SCIPgetBestSol(subscip);

         for( v = 0; v < nvars; ++v )
         {
            SCIP_VAR* subvar;
            SCIP_VAR* var;
            SCIP_Real fixval;
            SCIP_Bool infeasible;
            SCIP_Bool fixed;

            var = vars[v];

            subvar = (SCIP_VAR*)SCIPhashmapGetImage(varmapfw, var);
            fixval =  SCIPgetSolVal(subscip, sol, subvar);

            SCIP_CALL( SCIPfixVar(scip, var, fixval, &infeasible, &fixed) );
            assert(!infeasible);

            if( fixed )
               (*nfixedvars)++;
         }

         SCIP_CALL( SCIPdelCons(scip, cons) );
         (*ndelconss)++;

         break;
      }
      case SCIP_STATUS_NODELIMIT:
      case SCIP_STATUS_TIMELIMIT:
      case SCIP_STATUS_MEMLIMIT:
      {
         SCIP_VAR* var;
         SCIP_Bool infeasible;
         SCIP_Bool tightened;

         /* transfer the bound changes */
         for( v = 0; v < nvars; ++v )
         {
            var = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[v]);

            SCIP_CALL( SCIPtightenVarLb(scip, vars[v], SCIPvarGetLbGlobal(var), TRUE, &infeasible, &tightened) );
            assert(!infeasible);

            if( tightened )
               (*nchgbds)++;

            SCIP_CALL( SCIPtightenVarUb(scip, vars[v], SCIPvarGetUbGlobal(var), TRUE, &infeasible, &tightened) );
            assert(!infeasible);

            if( tightened )
               (*nchgbds)++;

         }
         break;
      }
      case SCIP_STATUS_UNKNOWN:
      case SCIP_STATUS_USERINTERRUPT:
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_GAPLIMIT:
      case SCIP_STATUS_SOLLIMIT:
      case SCIP_STATUS_BESTSOLLIMIT:
         SCIPerrorMessage("invalid status code <%d>\n", SCIPgetStatus(scip));
         return SCIP_INVALIDDATA;
      }
   }

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}

/** check if the capacity requirements can be satisfied */
static
SCIP_RETCODE checkCapacityRequirements(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is infeasible */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *cutoff = FALSE;
   v = 0;

   while( v < consdata->nvars )
   {
      if( consdata->demands[v] > consdata->capacity )
      {
         *cutoff = TRUE;
         break;
      }
      else if( consdata->demands[v] == 0 || consdata->durations[v] == 0 )
      {
         SCIP_CALL( deletePos(scip, cons, v) );
      }
      else
         v++;
   }

   return SCIP_OKAY;
}

/** propagate the cumulative condition */
static
SCIP_RETCODE propagateCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_CONS*            cons,               /**< constraint which is propagated (needed to SCIPinferVar**Cons()) */
   SCIP_Bool             usebinvars,         /**< is the binary representation used --> holes can be propagated */
   SCIP_Bool             usecoretimes,       /**< should core times be propagated */
   SCIP_Bool             usecoretimesholes,  /**< should core times be propagated to detect holes? */
   SCIP_Bool             useedgefinding,     /**< should edge finding be performed */
   SCIP_Bool             useenergeticreasoning,     /**< should energetic reasoning be performed */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   SCIP_Bool*            redundant,          /**< pointer to store if the constraint is redundant */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is infeasible */
   )
{
   assert(nchgbds != NULL);
   assert(initialized != NULL);
   assert(cutoff != NULL);
   assert((*cutoff) == FALSE);

   /**@todo avoid always sorting the variable array */

   /* check if the constraint is redundant */
   SCIP_CALL( consCheckRedundancy(scip, nvars, vars, durations, demands, capacity, redundant) );

   if( *redundant )
      return SCIP_OKAY;

   /* propagate the job cores until nothing else can be detected */
   if( usecoretimes )
   {
      SCIP_CALL( propagateCores(scip, nvars, vars, durations, demands, capacity, cons,
            nchgbds, initialized, cutoff) );

      if( *cutoff )
         return SCIP_OKAY;
   }

   /* check whether propagating the cores and creating holes helps */
   if( usebinvars && usecoretimesholes )
   {
      /* experimentally inefficient, but possible to be turned on */
      SCIP_CALL( propagateCoresForHoles(scip, nvars, vars, durations, demands, capacity, cons,
            nchgbds, initialized, cutoff) );

      if( *cutoff )
         return SCIP_OKAY;
   }

   if( useedgefinding )
   {
      /* check for overload, which may result in a cutoff */
      SCIP_CALL( checkOverload(scip, nvars, vars, durations, demands, capacity, cons,
            initialized, cutoff) );

      if( *cutoff )
         return SCIP_OKAY;

      /* perform edge-finding for lower bounds */
      SCIP_CALL( performEdgeFindingDetection(scip, TRUE, nvars, vars, durations, demands, capacity, cons,
            nchgbds, initialized, cutoff) );

      if( *cutoff )
         return SCIP_OKAY;

      /* perform edge-finding for upper bounds */
      SCIP_CALL( performEdgeFindingDetection(scip, FALSE, nvars, vars, durations, demands, capacity, cons,
            nchgbds, initialized, cutoff) );

      if( *cutoff )
         return SCIP_OKAY;
   }

   if( useenergeticreasoning )
   {
      /* perform energetic reasoning */
      SCIP_CALL( performEnergeticReasoning(scip, nvars, vars, durations, demands, capacity, cons,
            nchgbds, initialized, cutoff) );
   }

   return SCIP_OKAY;
}

/** propagate the cumulative constraint */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to propagate */
   SCIP_Bool             usebinvars,         /**< is the binary representation used --> holes can be propagated */
   SCIP_Bool             usecoretimes,       /**< should core times be propagated */
   SCIP_Bool             usecoretimesholes,  /**< should core times be propagated to detect holes? */
   SCIP_Bool             useedgefinding,     /**< should edge finding be performed */
   SCIP_Bool             useenergeticreasoning,     /**< should energetic reasoning be performed? */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   int*                  ndelconss,          /**< pointer to store the number of deleted constraints */
   SCIP_Bool*            cutoff              /**< pointer to store if the constraint is infeasible */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool initialized;
   SCIP_Bool redundant;
   int oldnchgbds;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   oldnchgbds = *nchgbds;
   initialized = FALSE;
   redundant = FALSE;

   SCIP_CALL( propagateCumulativeCondition(scip,
         consdata->nvars, consdata->vars, consdata->durations, consdata->demands, consdata->capacity, cons,
         usebinvars, usecoretimes, usecoretimesholes, useedgefinding, useenergeticreasoning,
         nchgbds, &redundant, &initialized, cutoff) );

   if( redundant )
   {
      SCIPdebugMessage("%s deletes cumulative constraint <%s> since it is redundant\n",
         SCIPgetDepth(scip) == 0 ? "globally" : "locally", SCIPconsGetName(cons));

      SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      (*ndelconss)++;
   }
   else
   {
      if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && initialized )
      {
         /* run conflict analysis since it was initialized */
         assert(*cutoff == TRUE);
         SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
      }

      /* if successful, reset age of constraint */
      if( *cutoff || *nchgbds > oldnchgbds )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }

   return SCIP_OKAY;
}

/** checks constraint for violation, and adds it as a cut if possible */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint to be separated */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   SCIP_Bool*            reducedom,          /**< pointer to store TRUE, if a domain was reduced */
   SCIP_Bool*            separated           /**< pointer to store TRUE, if a cut was found */
   )
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_ROW* row;
   SCIP_Real minfeasibility;
   SCIP_Bool useall;
   int ncuts;
   int r;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("separate cumulative constraint <%s>\n", SCIPconsGetName(cons));

   if( consdata->demandrows == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons, FALSE) );
   }

   minfeasibility = SCIPinfinity(scip);
   row = NULL;
   useall = FALSE;
   ncuts = 0;

   /* check each row that is not contained in LP */
   for( r = 0; r < consdata->ndemandrows; ++r )
   {
      if( !SCIProwIsInLP(consdata->demandrows[r]) )
      {
         SCIP_Real feasibility;

         if( sol != NULL )
            feasibility = SCIPgetRowSolFeasibility(scip, consdata->demandrows[r], sol);
         else
            feasibility = SCIPgetRowLPFeasibility(scip, consdata->demandrows[r]);

         if( useall )
         {
            if( SCIPisFeasNegative(scip, feasibility) )
            {
               SCIP_CALL( SCIPaddCut(scip, sol,  consdata->demandrows[r], FALSE) );
               ncuts++;
            }
         }
         else
         {
            if( minfeasibility > feasibility )
            {
               minfeasibility = feasibility;
               row = consdata->demandrows[r];
            }
         }
      }
   }

   if( !useall && SCIPisFeasNegative(scip, minfeasibility)  )
   {
      SCIPdebugMessage("cumulative constraint <%s> separated cut with feasibility <%g>\n",
         SCIPconsGetName(cons), minfeasibility);

      assert(row != NULL);
      SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );

      /* if successful, reset age of constraint */
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      (*separated) = TRUE;
   }
   else if( ncuts > 0 )
   {
      /* if successful, reset age of constraint */
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      (*separated) = TRUE;
   }

   return SCIP_OKAY;
}

/** checks constraint for violation, and adds it as a cut if possible */
static
SCIP_RETCODE separateCoverCutsCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logic or constraint to be separated */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool*            separated           /**< pointer to store TRUE, if a cut was found */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_ROW* row;
   SCIP_Real minfeasibility;
   int r;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("separate cumulative constraint <%s>\n", SCIPconsGetName(cons));

   /* collect the linking constraints */
   if( consdata->linkingconss == NULL )
   {
      SCIP_CALL( collectLinkingCons(scip, consdata) );
   }

   if( !consdata->covercuts )
   {
      SCIP_CALL( createCoverCuts(scip, cons) );
   }

   row = NULL;
   minfeasibility = SCIPinfinity(scip);

   /* check each row of small covers that is not contained in LP */
   for( r = 0; r < consdata->nscoverrows; ++r )
   {
      if( !SCIProwIsInLP(consdata->scoverrows[r]) )
      {
         SCIP_Real feasibility;

         assert(consdata->scoverrows[r] != NULL);
         if( sol != NULL )
            feasibility = SCIPgetRowSolFeasibility(scip, consdata->scoverrows[r], sol);
         else
            feasibility = SCIPgetRowLPFeasibility(scip, consdata->scoverrows[r]);

         if( minfeasibility > feasibility )
         {
            minfeasibility = feasibility;
            row =  consdata->scoverrows[r];
         }
      }
   }

   if( SCIPisFeasNegative(scip, minfeasibility) )
   {
      SCIPdebugMessage("cumulative constraint <%s> separated 1 cover cut with feasibility %g\n",
         SCIPconsGetName(cons), minfeasibility);

      assert(row != NULL);
      SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );

      /* if successful, reset age of constraint */
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      (*separated) = TRUE;
   }

   minfeasibility = SCIPinfinity(scip);
   row = NULL;

   /* check each row of small covers that is not contained in LP */
   for( r = 0; r < consdata->nbcoverrows; ++r )
   {
      if( !SCIProwIsInLP(consdata->bcoverrows[r]) )
      {
         SCIP_Real feasibility;

         assert(consdata->bcoverrows[r] != NULL);
         if( sol != NULL )
            feasibility = SCIPgetRowSolFeasibility(scip, consdata->bcoverrows[r], sol);
         else
            feasibility = SCIPgetRowLPFeasibility(scip, consdata->bcoverrows[r]);

         if( minfeasibility > feasibility )
         {
            minfeasibility = feasibility;
            row =  consdata->bcoverrows[r];
         }
      }
   }

   if( SCIPisFeasNegative(scip, minfeasibility) )
   {
      SCIPdebugMessage("cumulative constraint <%s> separated 1 cover cut with feasibility %g\n",
         SCIPconsGetName(cons), minfeasibility);

      assert(row != NULL);
      SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );

      /* if successful, reset age of constraint */
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      (*separated) = TRUE;
   }

   return SCIP_OKAY;
}

/** collect all integer variable which belong to jobs which can run at the point of interest */
static
SCIP_RETCODE collectIntVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_VAR***           activevars,         /**< jobs that are currently running */
   int*                  startindices,       /**< permutation with respect to the start times */
   int                   curtime,            /**< current point in time */
   int                   nstarted,           /**< number of jobs that start before the curtime or at curtime */
   int                   nfinished,          /**< number of jobs that finished before curtime or at curtime */
   SCIP_Bool             lower,              /**< shall cuts be created due to lower or upper bounds? */
   int*                  lhs                 /**< lhs for the new row sum of lbs + minoffset */
   )
{
   SCIP_VAR* var;
   int startindex;
   int endtime;
   int duration;
#ifndef NDEBUG
   int demand;
#endif
   int starttime;

   int varidx;
   int sumofstarts;
   int mindelta;
   int counter;

   counter = 0;
   sumofstarts = 0;

   mindelta = INT_MAX;

   startindex = nstarted - 1;

   /* search for the (nstarted - nfinished) jobs which are active at curtime */
   while( nstarted - nfinished > counter )
   {
      assert(startindex >= 0);

      /* collect job information */
      varidx = startindices[startindex];
      assert(varidx >= 0 && varidx < consdata->nvars);

      var = consdata->vars[varidx];
      duration = consdata->durations[varidx];
      assert(duration > 0);
#ifndef NDEBUG
      demand = consdata->demands[varidx];
      assert(demand > 0);
#endif
      assert(var != NULL);

      starttime = lower ? convertBoundToInt(scip, SCIPvarGetLbLocal(var)) : convertBoundToInt(scip, SCIPvarGetUbLocal(var));
      endtime = starttime + duration;

      /* check the end time of this job is larger than the curtime; in this case the job is still running */
      if( endtime > curtime )
      {
         (*activevars)[counter] = var;
         sumofstarts += starttime;
         mindelta = MIN(mindelta, endtime - curtime);
         counter++;
      }

      startindex--;
   }

   assert(mindelta > 0);
   *lhs = lower ? sumofstarts + mindelta : sumofstarts - mindelta;

   return SCIP_OKAY;
}

/** initialize the sorted event point arrays */
static
void createSortedEventpointsSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   int*                  starttimes,         /**< array to store sorted start events */
   int*                  endtimes,           /**< array to store sorted end events */
   int*                  startindices,       /**< permutation with respect to the start times */
   int*                  endindices,         /**< permutation with respect to the end times */
   int*                  nvars,              /**< number of variables that are integral */
   SCIP_Bool             lower               /**< shall the constraints be derived for lower or upper bounds? */
   )
{
   SCIP_VAR* var;
   int tmpnvars;
   int j;

   tmpnvars = consdata->nvars;
   *nvars = 0;

   /* assign variables, start and endpoints to arrays */
   for( j = 0; j < tmpnvars; ++j )
   {
      var = consdata->vars[j];
      assert(var != NULL);

      if( lower )
      {
         /* only consider jobs that are at their lower or upper bound */
         if( !SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, var))
            || !SCIPisFeasEQ(scip, SCIPgetSolVal(scip, sol, var), SCIPvarGetLbLocal(var)) )
            continue;

         if( consdata->durations[j] == 0 || consdata->demands[j] == 0 )
            continue;

         starttimes[*nvars] = convertBoundToInt(scip, SCIPgetSolVal(scip, sol, var));
         startindices[*nvars] = j;

         endtimes[*nvars] =  starttimes[*nvars] + consdata->durations[j];
         endindices[*nvars] = j;

         (*nvars) = *nvars + 1;

         SCIPdebugMessage("lower bounds are considered:\n");
         SCIPdebugMessage("%d: job[%d] starttime %d, endtime = %d, demand = %d\n ", *nvars-1,
            startindices[*nvars-1], starttimes[*nvars-1], starttimes[*nvars-1] + consdata->durations[startindices[*nvars-1]],
            consdata->demands[startindices[*nvars-1]]);
      }
      else
      {
         if( !SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, var))
            || !SCIPisFeasEQ(scip, SCIPgetSolVal(scip, sol, var), SCIPvarGetUbLocal(var)) )
            continue;

         starttimes[*nvars] = convertBoundToInt(scip, SCIPgetSolVal(scip, sol, var));
         startindices[*nvars] = j;

         endtimes[*nvars] =  starttimes[*nvars] + consdata->durations[j];
         endindices[*nvars] = j;

         (*nvars) = *nvars + 1;

         SCIPdebugMessage("upper bounds are considered:\n");
         SCIPdebugMessage("%d: job[%d] starttime %d, endtime = %d, demand = %d\n ", *nvars-1,
            startindices[*nvars-1], starttimes[*nvars-1], starttimes[*nvars-1] + consdata->durations[startindices[*nvars-1]],
            consdata->demands[startindices[*nvars-1]]);
      }
   }

   /* sort the arrays not-decreasing according to startsolvalues and endsolvalues (and sort the indices in the same way) */
   SCIPsortIntInt(starttimes, startindices, *nvars);
   SCIPsortIntInt(endtimes, endindices, *nvars);

#ifdef SCIP_DEBUG
   SCIPdebugMessage("sorted output\n");
   for( j = 0; j < *nvars; ++j )
   {
      SCIPdebugMessage("%d: job[%d] starttime %d, endtime = %d, demand = %d\n", j,
         startindices[j], starttimes[j], starttimes[j] + consdata->durations[startindices[j]],
         consdata->demands[startindices[j]]);
   }

   for( j = 0; j < *nvars; ++j )
   {
      SCIPdebugMessage("%d: job[%d] endtime %d,  demand = %d\n", j, endindices[j], endtimes[j],
         consdata->demands[endindices[j]]);
   }
   SCIPdebugMessage("capacity = %d\n", consdata->capacity);
#endif
}


/** this method creates a row for time point curtime which ensures the capacity restriction of the cumulative constraint */
static
SCIP_RETCODE createCapacityRestrictionIntvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be checked */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   int*                  startindices,       /**< permutation with respect to the start times */
   int                   curtime,            /**< current point in time */
   int                   nstarted,           /**< number of jobs that start before the curtime or at curtime */
   int                   nfinished,          /**< number of jobs that finished before curtime or at curtime */
   SCIP_Bool             lower               /**< shall cuts be created due to lower or upper bounds? */
   )
{
   SCIP_CONSDATA* consdata;
   char name[SCIP_MAXSTRLEN];
#ifndef NDEBUG
   int capacity;
#endif
   int lhs; /* left hand side of constraint */

   SCIP_VAR** activevars;
   SCIP_ROW* row;

   int v;

   assert(nstarted > nfinished);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars > 0);

#ifndef NDEBUG
   capacity = consdata->capacity;
   assert(capacity > 0);
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &activevars, nstarted-nfinished) );

   SCIP_CALL( collectIntVars(scip, consdata, &activevars,
         startindices, curtime, nstarted, nfinished, lower, &lhs ) );

   if( lower )
   {
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "lower(%d)", curtime);

      SCIP_CALL( SCIPcreateEmptyRow(scip, &row, name, (SCIP_Real) lhs, SCIPinfinity(scip),  TRUE, FALSE, SCIPconsIsRemovable(cons)) );
   }
   else
   {
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "upper(%d)", curtime);
      SCIP_CALL( SCIPcreateEmptyRow(scip, &row, name, -SCIPinfinity(scip), (SCIP_Real) lhs, TRUE, FALSE, SCIPconsIsRemovable(cons)) );
   }
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   for( v = 0; v < nstarted - nfinished; ++v )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, row, activevars[v], 1.) );
   }

   SCIP_CALL( SCIPflushRowExtensions(scip, row) );
   SCIPdebug( SCIP_CALL(SCIPprintRow(scip, row, NULL)) );

   SCIP_CALL( SCIPaddCut(scip, sol, row, TRUE) );

   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   /* free buffers */
   SCIPfreeBufferArrayNull(scip, &activevars);

   return SCIP_OKAY;
}

/** checks constraint for violation, and adds it as a cut if possible */
static
SCIP_RETCODE separateConsOnIntegerVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< cumulative constraint to be separated */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool             lower,              /**< shall cuts be created according to lower bounds? */
   SCIP_Bool*            separated           /**< pointer to store TRUE, if a cut was found */
   )
{

   SCIP_CONSDATA* consdata;

   int* starttimes;         /* stores when each job is starting */
   int* endtimes;           /* stores when each job ends */
   int* startindices;       /* we will sort the startsolvalues, thus we need to know which index of a job it corresponds to */
   int* endindices;         /* we will sort the endsolvalues, thus we need to know which index of a job it corresponds to */

   int nvars;               /* number of activities for this constraint */
   int freecapacity;        /* remaining capacity */
   int curtime;             /* point in time which we are just checking */
   int endindex;            /* index of endsolvalues with: endsolvalues[endindex] > curtime */

   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   /* if no activities are associated with this cumulative then this constraint is redundant */
   if( nvars == 0 )
      return SCIP_OKAY;

   assert(consdata->vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &starttimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endtimes, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &startindices, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &endindices, nvars) );

   SCIPdebugMessage("create sorted event points for cumulative constraint <%s> with %d jobs\n",
      SCIPconsGetName(cons), nvars);

   /* create event point arrays */
   createSortedEventpointsSol(scip, consdata, sol, starttimes, endtimes, startindices, endindices, &nvars, lower);

   endindex = 0;
   freecapacity = consdata->capacity;

   /* check each startpoint of a job whether the capacity is kept or not */
   /* only check those 'nvars' that are not fractional (only those were sorted in 'createsortedeventpointssol') */
   for( j = 0; j < nvars; ++j )
   {
      curtime = starttimes[j];

      /* remove the capacity requirements for all job which start at the curtime */
      subtractStartingJobDemands(consdata, curtime, starttimes, startindices, &freecapacity, &j, nvars);

      /* add the capacity requirements for all job which end at the curtime */
      addEndingJobDemands(consdata, curtime, endtimes, endindices, &freecapacity, &endindex, nvars);

      assert(freecapacity <= consdata->capacity);
      assert(endindex <= nvars);

      /* endindex - points to the next job which will finish */
      /* j - points to the last job that has been released */

      /* if free capacity is smaller than zero, then add rows to the LP */
      if( freecapacity < 0 )
      {
         /* create capacity restriction row for current event point */
         SCIP_CALL( createCapacityRestrictionIntvars(scip, cons, sol, startindices, curtime, j+1, endindex, lower) );
         *separated = TRUE;
      }
   } /*lint --e{850}*/

   /* free all buffer arrays */
   SCIPfreeBufferArray(scip, &endindices);
   SCIPfreeBufferArray(scip, &startindices);
   SCIPfreeBufferArray(scip, &endtimes);
   SCIPfreeBufferArray(scip, &starttimes);

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyCumulative)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrCumulative(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeCumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdataFree(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** initialization method of constraint handler (called after problem was transformed) */
#define consInitCumulative NULL

/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitCumulative NULL

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreCumulative)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   (*result) = SCIP_FEASIBLE;

   /* check all constraints for trivial feasibility  or infeasibility */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert(cons != NULL);

      if( !checkDemands(scip, cons) )
      {
         (*result) = SCIP_CUTOFF;
         break;
      }
   }

   return SCIP_OKAY;
}

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreCumulative NULL

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolCumulative NULL

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolCumulative)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* free rows */
      SCIP_CALL( consdataFreeRows(scip, &consdata) );
   }

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteCumulative)
{  /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL );
   assert(*consdata != NULL );

   /* free cumulative constraint data */
   SCIP_CALL( consdataFree(scip, consdata ) );

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransCumulative)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   assert(conshdlr != NULL);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->demandrows == NULL);

   SCIPdebugMessage("transform cumulative constraint <%s>\n", SCIPconsGetName(sourcecons));

   /* create constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->vars, sourcedata->linkingconss,
         sourcedata->durations, sourcedata->demands, sourcedata->nvars, sourcedata->capacity) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpCumulative)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdebugMessage("initialize LP relaxation for %d cumulative constraints\n", nconss);

   if( conshdlrdata->usebinvars )
   {
      /* add rows to LP */
      for( c = 0; c < nconss; ++c )
      {
         assert(SCIPconsIsInitial(conss[c]));
         SCIP_CALL( addRelaxation(scip, conss[c], conshdlrdata->cutsasconss) );

         if( conshdlrdata->cutsasconss )
         {
            SCIP_CALL( SCIPrestartSolve(scip) );
         }
      }
   }

   /**@todo if we want to use only the integer variables; only these will be in cuts
    *       create some initial cuts */

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpCumulative)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool reducedom;
   SCIP_Bool separated;
   int c;

   SCIPdebugMessage("consSepalpCumulative\n");

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdebugMessage("separating %d/%d cumulative constraints\n", nusefulconss, nconss);

   cutoff = FALSE;
   reducedom = FALSE;
   separated = FALSE;
   (*result) = SCIP_DIDNOTFIND;

   if( conshdlrdata->usebinvars )
   {
      /* check all useful cumulative constraints for feasibility  */
      for( c = 0; c < nusefulconss && !reducedom && !cutoff; ++c )
      {
         SCIP_CALL( separateCons(scip, conss[c], NULL, &cutoff, &reducedom, &separated) );
      }

      if( !cutoff && !reducedom && conshdlrdata->usecovercuts )
      {
         for( c = 0; c < nusefulconss; ++c )
         {
            SCIP_CALL( separateCoverCutsCons(scip, conss[c], NULL, &separated) );
         }
      }
   }
   else
   {
      /* separate cuts containing only integer variables */
      for( c = 0; c < nusefulconss; ++c )
      {
         SCIP_CALL( separateConsOnIntegerVariables(scip, conss[c], NULL, TRUE, &separated) );
         SCIP_CALL( separateConsOnIntegerVariables(scip, conss[c], NULL, FALSE, &separated) );
      }
   }

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reducedom )
      *result = SCIP_REDUCEDDOM;
   else if( separated )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolCumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool reducedom;
   SCIP_Bool separated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   SCIPdebugMessage("separating %d/%d cumulative constraints\n", nusefulconss, nconss);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   reducedom = FALSE;
   separated = FALSE;
   (*result) = SCIP_DIDNOTFIND;

   if( conshdlrdata->usebinvars )
   {
      /* check all useful cumulative constraints for feasibility  */
      for( c = 0; c < nusefulconss && !cutoff && !reducedom; ++c )
      {
         SCIP_CALL( separateCons(scip, conss[c], NULL, &cutoff, &reducedom, &separated) );
      }

      if( !cutoff && !reducedom && conshdlrdata->usecovercuts )
      {
         for( c = 0; c < nusefulconss; ++c )
         {
            SCIP_CALL( separateCoverCutsCons(scip, conss[c], sol, &separated) );
         }
      }
   }
   else
   {
      /* separate cuts containing only integer variables */
      for( c = 0; c < nusefulconss; ++c )
      {
         SCIP_CALL( separateConsOnIntegerVariables(scip, conss[c], sol, TRUE, &separated) );
         SCIP_CALL( separateConsOnIntegerVariables(scip, conss[c], sol, FALSE, &separated) );
      }
   }

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( reducedom )
      *result = SCIP_REDUCEDDOM;
   else if( separated )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpCumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   SCIP_Bool reducedom;
   SCIP_Bool separated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( solinfeasible )
   {
      *result = SCIP_INFEASIBLE;
      return SCIP_OKAY;
   }

   cutoff = FALSE;
   reducedom = FALSE;
   separated = FALSE;

   SCIPdebugMessage("LP enforcing %d useful resource constraints of %d constraints\n", nusefulconss, nconss);

   if( conshdlrdata->usebinvars )
   {

      /* check all useful cumulative constraints for feasibility */
      for( c = 0; c < nusefulconss && !cutoff && !reducedom; ++c )
      {
         SCIP_CALL( separateCons(scip, conss[c], NULL, &cutoff, &reducedom, &separated) );
      }

      /* check all obsolete cumulative constraints for feasibility */
      for( c = nusefulconss; c < nconss && !cutoff && !reducedom && !separated; ++c )
      {
         SCIP_CALL( separateCons(scip, conss[c], NULL, &cutoff, &reducedom, &separated) );
      }

      if( cutoff )
         *result = SCIP_CUTOFF;
      else if( reducedom )
         *result = SCIP_REDUCEDDOM;
      else if( separated )
         *result = SCIP_SEPARATED;
      else
         (*result) = SCIP_FEASIBLE;
   }
   else
   {
      /* it is no longer clear how to forbid a solution by cuts on integer variables -> only check solution */
      SCIP_Bool violated;

      violated = FALSE;

      for( c = 0; c < nconss && !violated; ++c )
      {
         SCIP_CALL( checkCons(scip, conss[c], NULL, &violated, FALSE) );
      }

      if( violated )
         *result = SCIP_INFEASIBLE;
      else
         *result = SCIP_FEASIBLE;
   }

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsCumulative)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int c;

   SCIPdebugMessage("method: enforce pseudo solution\n");

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   if( objinfeasible )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   violated = FALSE;

   (*result) = SCIP_FEASIBLE;
   for( c = 0; c < nconss && !violated; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], NULL, &violated, FALSE) );
   }

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckCumulative)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   violated = FALSE;

   for( c = 0; c < nconss && !violated; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], sol, &violated, printreason) );
   }

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropCumulative)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   int nchgbds;
   int ndelconss;
   int c;

   SCIPdebugMessage("propagate cumulative constraints\n");

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(nconss == 0 || conss != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   nchgbds = 0;
   ndelconss = 0;
   cutoff = FALSE;
   (*result) = SCIP_DIDNOTRUN;

   /* propagate all useful constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( propagateCons(scip, conss[c],
            conshdlrdata->usebinvars, conshdlrdata->usecoretimes, conshdlrdata->usecoretimesholes,
            conshdlrdata->useedgefinding, conshdlrdata->useenergeticreasoning,
            &nchgbds, &ndelconss, &cutoff) );
   }

   if( !cutoff && nchgbds == 0 )
   {
      /* propagate all other constraints */
      for( c = nusefulconss; c < nconss && !cutoff; ++c )
      {
         SCIP_CALL( propagateCons(scip, conss[c],
               conshdlrdata->usebinvars, conshdlrdata->usecoretimes, conshdlrdata->usecoretimesholes,
               conshdlrdata->useedgefinding, conshdlrdata->useenergeticreasoning,
               &nchgbds, &ndelconss, &cutoff) );
      }
   }

   if( cutoff )
   {
      SCIPdebugMessage("detected infeasible\n");
      *result = SCIP_CUTOFF;
   }
   else if( nchgbds > 0 )
   {
      SCIPdebugMessage("delete (locally) %d constraints and changed %d variable bounds\n", ndelconss, nchgbds);
      *result = SCIP_REDUCEDDOM;
   }
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolCumulative)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_Bool cutoff;
   SCIP_Bool unbounded;
   int oldnchgbds;
   int oldndelconss;
   int oldnfixedvars;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("presolve cumulative constraints\n");

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTRUN;

   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnfixedvars = * nfixedvars;
   cutoff = FALSE;
   unbounded = FALSE;

   /* in the first round check if all demands are smaller or equal to the capacity or have a energy of zero */
   if( nrounds == 0 )
   {
      for( c = 0; c < nconss && !cutoff; ++c )
      {
         SCIP_CALL( checkCapacityRequirements(scip, conss[c], &cutoff) );
      }
   }

   /* process constraints */
   for( c = 0; c < nconss && !cutoff && !unbounded; ++c )
   {
      cons = conss[c];
      assert(cons != NULL);

      /* remove irrelevant jobs */
      SCIP_CALL( removeIrrelevantJobs(scip, cons, nchgbds) );

      /* divide demands by their greatest common divisor */
      normalizeDemands(cons, nchgcoefs, nchgsides);

      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, cons,
            conshdlrdata->usebinvars, conshdlrdata->usecoretimes, conshdlrdata->usecoretimesholes,
            conshdlrdata->useedgefinding, TRUE,
            nchgbds, ndelconss, &cutoff) );

      if( !SCIPconsIsDeleted(cons) )
      {
         SCIP_CALL( dualPresolving(scip, cons, conshdlrdata->maxnodes, nchgbds, nfixedvars, ndelconss, &cutoff, &unbounded) );
      }
   }

   SCIPdebugMessage("cutoff %u delete %d constraints and changed %d variable bounds\n",
      cutoff, *ndelconss - oldndelconss, *nchgbds - oldnchgbds);

   /* evaluate the presolving round */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( unbounded )
      *result = SCIP_UNBOUNDED;
   else if( *nchgbds > oldnchgbds || *nfixedvars > oldnfixedvars || *ndelconss > oldndelconss )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropCumulative)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(infervar != NULL);
   assert(bdchgidx != NULL);

   SCIPdebugMessage("resolve propagation for variable <%s> and cumulative constraint <%s> with rule %d\n",
      SCIPvarGetName(infervar), SCIPconsGetName(cons), inferInfoGetProprule(intToInferInfo(inferinfo)));

   /* process constraint */
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( respropCumulativeCondition(scip, consdata->nvars, consdata->vars,
         consdata->durations, consdata->demands, consdata->capacity,
         infervar, intToInferInfo(inferinfo), boundtype, bdchgidx, result) );

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockCumulative)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( v = 0; v < consdata->nvars; ++v )
   {
      /* the integer start variable should not get rounded in both direction  */
      assert(consdata->vars[v] != NULL);
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[v], nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   return SCIP_OKAY;
}

/** constraint activation notification method of constraint handler */
#define consActiveCumulative NULL

/** constraint deactivation notification method of constraint handler */
#define consDeactiveCumulative NULL

/** constraint enabling notification method of constraint handler */
#define consEnableCumulative NULL

/** constraint disabling notification method of constraint handler */
#define consDisableCumulative NULL

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintCumulative)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdataPrint(scip, SCIPconsGetData(cons), file);

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyCumulative)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourceconsdata;
   SCIP_VAR** sourcevars;
   SCIP_VAR** vars;
   const char* consname;

   int nvars;
   int v;

   sourceconsdata = SCIPconsGetData(sourcecons);
   assert(sourceconsdata != NULL);

   /* get variables of the source constraint */
   nvars = sourceconsdata->nvars;
   sourcevars = sourceconsdata->vars;

   (*valid) = TRUE;

   if( nvars == 0 )
      return SCIP_OKAY;

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

   for( v = 0; v < nvars && *valid; ++v )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[v], &vars[v], varmap, consmap, global, valid) );
      assert(!(*valid) || vars[v] != NULL);
   }

   /* only create the target constraint, if all variables could be copied */
   if( *valid )
   {
      if( name != NULL )
         consname = name;
      else
         consname = SCIPconsGetName(sourcecons);

      /* copy the logic using the linear constraint copy method */
      SCIP_CALL( SCIPcreateConsCumulative(scip, cons, consname, nvars, vars,
            sourceconsdata->durations, sourceconsdata->demands, sourceconsdata->capacity,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   /* free buffer array */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** variable deletion method of constraint handler */
#define consDelvarsCumulative NULL


/** constraint parsing method of constraint handler */
#define consParseCumulative NULL

/*
 * constraint specific interface methods
 */

/** creates the handler for cumulative constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrCumulative(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create cumulative constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING,
         conshdlrCopyCumulative,
         consFreeCumulative, consInitCumulative, consExitCumulative,
         consInitpreCumulative, consExitpreCumulative, consInitsolCumulative, consExitsolCumulative,
         consDeleteCumulative, consTransCumulative, consInitlpCumulative,
         consSepalpCumulative, consSepasolCumulative, consEnfolpCumulative, consEnfopsCumulative, consCheckCumulative,
         consPropCumulative, consPresolCumulative, consRespropCumulative, consLockCumulative,
         consActiveCumulative, consDeactiveCumulative,
         consEnableCumulative, consDisableCumulative, consDelvarsCumulative,
         consPrintCumulative, consCopyCumulative, consParseCumulative,
         conshdlrdata) );

   /* set default values for constraint handler data */
   conshdlrdata->lastsepanode = -1;

   /* add cumulative constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/usebinvars", "should the binary representation be used?",
         &conshdlrdata->usebinvars, FALSE, DEFAULT_USEBINVARS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/usecoretimes", "should coretimes be propagated?",
         &conshdlrdata->usecoretimes, FALSE, DEFAULT_USECORETIMES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/usecoretimesholes", "should coretimes be propagated to detect holes?",
         &conshdlrdata->usecoretimesholes, FALSE, DEFAULT_USECORETIMESHOLES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/localcuts", "should cuts be added only locally?",
         &conshdlrdata->localcuts, FALSE, DEFAULT_LOCALCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/usecovercuts", "should covering cuts be added every node?",
         &conshdlrdata->usecovercuts, FALSE, DEFAULT_USECOVERCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/useedgefinding", "should edge finding be used?",
         &conshdlrdata->useedgefinding, FALSE, DEFAULT_USEEDGEFINDING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/useenergeticreasoning", "should energetic reasoning be used?",
         &conshdlrdata->useenergeticreasoning, FALSE, DEFAULT_USEENERGETICREASONING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/cutsasconss",
         "should the cumulative constraint create cuts as knapsack constraints?",
         &conshdlrdata->cutsasconss, FALSE, DEFAULT_CUTSASCONSS, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip,
         "constraints/"CONSHDLR_NAME"/maxnodes",
         "maximum nodes to solve an independent cumulative constraint (-1: no limit)",
         &conshdlrdata->maxnodes, TRUE, DEFAULT_MAXNODES, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a cumulative constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);

   /* find the precedence constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage(""CONSHDLR_NAME" constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   SCIPdebugMessage("create cumulative constraint <%s> with %d jobs\n", name, nvars);

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, vars, NULL, durations, demands, nvars, capacity) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata,
         initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** returns the activities of the cumulative constraint */
SCIP_VAR** SCIPgetVarsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cumulative constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** returns the activities of the cumulative constraint */
int SCIPgetNVarsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cumulative constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** returns the capacity of the cumulative constraint */
int SCIPgetCapacityCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cumulative constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->capacity;
}

/** returns the durations of the cumulative constraint */
int* SCIPgetDurationsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cumulative constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->durations;
}

/** returns the demands of the cumulative constraint */
int* SCIPgetDemandsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cumulative constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->demands;
}

/** check for the given starting time variables with their demands and durations if the cumulative conditions for the
 *  given solution is satisfied
 */
SCIP_RETCODE SCIPcheckCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
   SCIP_Bool*            violated,           /**< pointer to store if the cumulative condition is violated */
   SCIP_CONS*            cons,               /**< constraint which is checked */
   SCIP_Bool             printreason         /**< should the reason for the violation be printed? */
   )
{
   assert(scip != NULL);
   assert(violated != NULL);

   SCIP_CALL( checkCumulativeCondition(scip, sol, nvars, vars, durations, demands, capacity,
         violated, cons, printreason) );

   return SCIP_OKAY;
}

/** propagate the given cumulative condition */
SCIP_RETCODE SCIPpropCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
   SCIP_CONS*            cons,               /**< constraint which gets propagated */
   int*                  nchgbds,            /**< pointer to store the number of variable bound changes */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            cutoff              /**< pointer to store if the cumulative condition is violated */
   )
{
   SCIP_Bool redundant;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(initialized != NULL);
   assert(*initialized == FALSE);
   assert(cutoff != NULL);
   assert(*cutoff == FALSE);

   redundant = FALSE;

   SCIP_CALL( propagateCumulativeCondition(scip,nvars, vars, durations, demands, capacity, cons,
         FALSE, TRUE, FALSE, TRUE, FALSE, nchgbds, &redundant, initialized, cutoff) );

   return SCIP_OKAY;
}

/** resolve propagation w.r.t. the cumulative condition */
SCIP_RETCODE SCIPrespropCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_VAR*             infervar,           /**< the conflict variable whose bound change has to be resolved */
   int                   inferinfo,          /**< the user information */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   )
{
   SCIP_CALL( respropCumulativeCondition(scip, nvars, vars, durations, demands, capacity,
         infervar, intToInferInfo(inferinfo), boundtype, bdchgidx, result) );

   return SCIP_OKAY;
}

/** create a new cumulative profile for the given capacity */
SCIP_RETCODE SCIPprofileCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   CUMULATIVEPROFILE**   profile,            /**< pointer to store the create profile */
   int                   capacity,           /**< capacity for this profile */
   int                   maxtimepoints       /**< maximum number of time points */
   )
{
   assert(scip != NULL);
   assert(profile != NULL);
   assert(capacity > 0);
   assert(maxtimepoints > 0);

   /* initialize memory */
   SCIP_CALL( SCIPallocMemory(scip, profile) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*profile)->timepoints, maxtimepoints) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*profile)->freecapacities, maxtimepoints) );

   /* set up cumulative profile for use */
   (*profile)->ntimepoints = 2;
   (*profile)->timepoints[0] = 0;
   (*profile)->timepoints[1] = INT_MAX;
   (*profile)->freecapacities[0] = capacity;
   (*profile)->freecapacities[1] = 0;
   (*profile)->arraysize = maxtimepoints;

   return SCIP_OKAY;
}

/** frees given profile */
void SCIPprofileFree(
   SCIP*                 scip,               /**< SCIP data structure */
   CUMULATIVEPROFILE**   profile             /**< pointer to the profile */
   )
{
   assert(scip != NULL);
   assert(profile != NULL);

   /* free memory */
   SCIPfreeMemoryArray(scip, &(*profile)->timepoints);
   SCIPfreeMemoryArray(scip, &(*profile)->freecapacities);
   SCIPfreeMemory(scip, profile);
}

/** resizes the cumulative profile array */
SCIP_RETCODE SCIPprofileResize(
   SCIP*                 scip,               /**< SCIP data structure */
   CUMULATIVEPROFILE*    profile,            /**< cumulative profile to resize */
   int                   newminsize          /**< minimum size to ensure */
   )
{
   assert(scip != NULL);
   assert(profile != NULL);
   assert(newminsize >= 0);
   assert(profile->timepoints != NULL);
   assert(profile->freecapacities != NULL);

   if( profile->ntimepoints >= newminsize )
      return SCIP_OKAY;

   /* Grow arrays of times and free capacity */
   SCIP_CALL( SCIPreallocMemoryArray(scip, &profile->timepoints, newminsize) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &profile->freecapacities, newminsize) );
   profile->arraysize = newminsize;

   return SCIP_OKAY;
}

/** from the given job, the core time is computed. If core is non-empty the cumulative profile will be updated otherwise
 *  nothing happens
 */
void SCIPprofileInsertCore(
   SCIP*                 scip,               /**< SCIP data structure */
   CUMULATIVEPROFILE*    profile,            /**< profile to use */
   SCIP_VAR*             var,                /**< integer variable which corresponds to the starting point of the job */
   int                   duration,           /**< duration of the job */
   int                   demand,             /**< demand of the job */
   SCIP_Bool*            core,               /**< pointer to store if the corresponds job has a core */
   SCIP_Bool*            fixed,              /**< pointer to store if the job is fixed due to its bounds */
   SCIP_Bool*            infeasible          /**< pointer to store if the job does not fit due to capacity */
   )
{
   int begin;
   int end;
   int lb;
   int ub;

   assert(core != NULL);
   assert(fixed != NULL);
   assert(infeasible != NULL);

   (*infeasible) = FALSE;
   (*fixed) = FALSE;
   (*core) = FALSE;

   lb = convertBoundToInt(scip, SCIPvarGetLbLocal(var));
   ub = convertBoundToInt(scip, SCIPvarGetUbLocal(var));

   if( ub - lb == 0 )
      (*fixed) = TRUE;

   begin = ub;
   end = lb + duration;

   /* check if a core exists */
   if( begin < end )
   {
      /* job has a nonempty core and will be inserted */
      (*core) = TRUE;

      /* insert core into the profile */
#ifdef PROFILE_DEBUG
      SCIPdebugMessage("before inserting: \n");
      SCIPprofilePrintOut(profile);
      SCIPdebugMessage("insert core from var <%s>[%d,%d]: [%d,%d] [%d]\n", SCIPvarGetName(var), lb, ub, begin, end, demand);
#endif

      SCIPprofileUpdate(profile, begin, end, demand, infeasible);

#ifdef PROFILE_DEBUG
      {
         int i;
         SCIPdebugMessage("after inserting: %u\n", *infeasible);
         SCIPprofilePrintOut(profile);

         for( i =0; i < profile->ntimepoints-1; ++i )
         {
            assert(profile->timepoints[i] < profile->timepoints[i+1]);
         }
      }
#endif
   }
}

/** subtracts the demand from the profile during core time of the job */
void SCIPprofileDeleteCore(
   SCIP*                 scip,               /**< SCIP data structure */
   CUMULATIVEPROFILE*    profile,            /**< profile to use */
   SCIP_VAR*             var,                /**< integer variable which corresponds to the starting point of the job */
   int                   duration,           /**< duration of the job */
   int                   demand,             /**< demand of the job */
   SCIP_Bool*            core                /**< pointer to store if the corresponds job has a core, or NULL */
   )
{
   int begin;
   int end;
   SCIP_Bool infeasible;

   begin = convertBoundToInt(scip, SCIPvarGetUbLocal(var));
   end = convertBoundToInt(scip, SCIPvarGetLbLocal(var)) + duration;

   if( begin >= end )
   {
      if(core != NULL)
         *core = FALSE;

      return;
   }

   if( core != NULL )
      *core = TRUE;

#ifndef NDEBUG
   {
      /* check if the begin and end time points of the core correspond to a time point in the profile; this should be
       * the case since we added the core before to the profile */
      int pos;
      assert(SCIPprofileFindLowerBound(profile, begin, &pos));
      assert(SCIPprofileFindLowerBound(profile, end, &pos));
   }
#endif

   /* remove the core of the job from the current profile */
#ifdef PROFILE_DEBUG
   SCIPdebugMessage("before deleting:\n");
   SCIPprofilePrintOut(profile);

   SCIPdebugMessage("delete core from var <%s>: [%d,%d] [%d]\n",
      SCIPvarGetName(var), begin, end, demand);
#endif

   SCIPprofileUpdate(profile, begin, end, -demand, &infeasible);

#ifdef PROFILE_DEBUG
   SCIPdebugMessage("after deleting: %u\n", infeasible);
   SCIPprofilePrintOut(profile);
#endif
   assert(!infeasible);
}


/** output of the given profile */
void SCIPprofilePrint(
   SCIP*                 scip,               /**< SCIP data structure */
   CUMULATIVEPROFILE*    profile,            /**< profile to output */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int t;

   for( t = 0; t < profile->ntimepoints; ++t )
   {
      SCIPinfoMessage(scip, file, "i: %d, tp: %d, fc: %d ;", t, profile->timepoints[t], profile-> freecapacities[t]);
   }

   SCIPinfoMessage(scip, file,"\n");
}


/** return if the given time point exists in the profile and stores the position of the given time point if it exists;
 *  otherwise the position of the next smaller existing time point */
SCIP_Bool SCIPprofileFindLowerBound(
   CUMULATIVEPROFILE*    profile,              /**< profile to search in */
   int                   timepoint,            /**< time point to search for */
   int*                  pos                   /**< pointer to store the position */
   )
{
   assert(profile != NULL);
   assert(timepoint >= 0);
   assert(profile->ntimepoints > 0);
   assert(profile->timepoints[0] == 0);

   /* find the position of timepoint in the timepoints array via binary search */
   if( SCIPsortedvecFindInt(profile->timepoints, timepoint, profile->ntimepoints, pos) )
      return TRUE;

   assert(*pos > 0);
   (*pos)--;

   return FALSE;
}

/** inserts the given time point into the profile if it this time point does not exists yet; returns its position in the
 *  time point array */
int SCIPprofileInsertTimepoint(
   CUMULATIVEPROFILE*    profile,            /**< profile to insert the time point */
   int                   timepoint           /**< time point to insert */
   )
{
   int pos;
#ifndef NDEBUG
   int i;
#endif

   assert(profile != NULL);
   assert(timepoint >= 0);
   assert(profile->arraysize >= profile->ntimepoints);

   if( timepoint == 0 )
      return 0;

   /* get the position of the given time point in the profile array if it exists; otherwise the position of the next
    * smaller existing time point */
   if( SCIPprofileFindLowerBound(profile, timepoint, &pos) )
   {
      /* if the time point exists return the corresponding position */
      assert(pos >= 0 && pos < profile->ntimepoints);
      return pos;
   }

   assert(pos >= 0 && pos < profile->ntimepoints);
   assert(timepoint >= profile->timepoints[pos]);
   assert(pos + 1 < profile->arraysize);

   /* insert new time point into the (sorted) profile */
   SCIPsortedvecInsertIntInt(profile->timepoints, profile->freecapacities, timepoint, profile->freecapacities[pos], 
      &profile->ntimepoints, NULL);

#ifndef NDEBUG
   /* check if the time points are sorted */
   for( i = 1; i < profile->ntimepoints; ++i )
      assert(profile->timepoints[i-1] < profile->timepoints[i]);
#endif

   return pos+1;
}

/** updates the profile due to inserting and removing a new job */
void SCIPprofileUpdate(
   CUMULATIVEPROFILE*    profile,            /**< profile to update */
   int                   starttime,          /**< time point to start */
   int                   endtime,            /**< time point to end */
   int                   demand,             /**< demand of the job */
   SCIP_Bool*            infeasible          /**< pointer to store if the update is infeasible */
   )
{
   int startpos;
   int endpos;

   assert(profile != NULL);
   assert(infeasible != NULL);
   assert(profile->arraysize >= profile->ntimepoints);
   assert(starttime >= 0 && endtime >= starttime);

   (*infeasible) = FALSE;

   if( starttime == endtime )
      return;

   /* get position of the starttime in profile */
   startpos = SCIPprofileInsertTimepoint(profile, starttime);
   assert(profile->timepoints[startpos] == starttime);

   /* get position of the endtime in profile */
   endpos = SCIPprofileInsertTimepoint(profile, endtime);
   assert(profile->timepoints[endpos] == endtime );

   assert(startpos < endpos);
   assert(profile->arraysize >= profile->ntimepoints);

   /* remove/add the given demand from the profile */
   for( ; startpos < endpos; ++startpos )
   {
      profile->freecapacities[startpos] -= demand;

      if( profile->freecapacities[startpos] < 0 )
      {
         *infeasible = TRUE;
         break;
      }
   }
}

/** returns TRUE if the job (given by its  demand and duration) can be inserted at the given time point; otherwise FALSE */
SCIP_Bool SCIPprofileIsFeasibleStart(
   CUMULATIVEPROFILE*    profile,            /**< Cumulative profile to use */
   int                   timepoint,          /**< time point to start */
   int                   duration,           /**< duration of the job */
   int                   demand,             /**< the demand of the job */
   int*                  pos                 /**< pointer to store the earliest position where the job does not fit */
   )
{
   int endtime;
   int startpos;
   int endpos;
   int p;

   assert(profile != NULL);
   assert(timepoint >= 0);
   assert(demand >= 0);
   assert(pos != NULL);

   if( duration == 0 )
      return TRUE;

   endtime = timepoint + duration;

   /* check if the activity fits at timepoint */
   (void)SCIPprofileFindLowerBound(profile, timepoint, &startpos);

   if( !SCIPprofileFindLowerBound(profile, endtime, &endpos) )
      endpos++;

   assert(profile->timepoints[startpos] <= timepoint);
   assert(profile->timepoints[endpos] >= endtime);

   for( p = startpos; p < endpos; ++p )
   {
      if( profile->freecapacities[p] < demand )
      {
         (*pos) = p;
         return FALSE;
      }
   }

   return TRUE;
}

/** return the earliest possible starting point within the time interval [lb,ub] for a given job (given by its duration
 *  and demand) */
int SCIPprofileGetEarliestFeasibleStart(
   CUMULATIVEPROFILE*    profile,            /**< profile to use */
   int                   lb,                 /**< earliest possible start point */
   int                   ub,                 /**< latest possible start point */
   int                   duration,           /**< duration of the job */
   int                   demand,             /**< demand of the job */
   SCIP_Bool*            infeasible          /**< pointer store if the job cannot be scheduled */
   )
{
   int starttime;
   int pos;

   assert(profile != NULL);
   assert(lb >= 0);
   assert(duration >= 0);
   assert(demand >= 0);
   assert(infeasible != NULL);
   assert(profile->timepoints[profile->ntimepoints-1] > ub);

   if( lb > ub )
   {
      *infeasible = TRUE;
      return lb;
   }

   if( duration == 0 || demand == 0 )
   {
      *infeasible = FALSE;
      return lb;
   }

   starttime = lb;

   (void)SCIPprofileFindLowerBound(profile, starttime, &pos);
   assert(profile->timepoints[pos] <= starttime);

   (*infeasible) = TRUE;

   while( (*infeasible) && starttime <= ub )
   {
      if( SCIPprofileIsFeasibleStart(profile, starttime, duration, demand, &pos) )
      {
         (*infeasible) = FALSE;
         return starttime;
      }

      /* the job did not fit into the profile since at time point "pos" not enough capacity is available; therefore we
       * can proceed with the next time point  */
      assert(profile->freecapacities[pos] < demand);
      pos++;

      /* check if we exceed the time point array */
      if( pos >= profile->ntimepoints )
         break;

      starttime = profile->timepoints[pos];
   }

   assert(*infeasible || starttime <= ub);
   return starttime;
}

/** return the latest possible starting point within the time interval [lb,ub] for a given job (given by its duration
 *  and demand) */
int SCIPprofileGetLatestFeasibleStart(
   CUMULATIVEPROFILE*    profile,            /**< profile to use */
   int                   lb,                 /**< earliest possible start point */
   int                   ub,                 /**< latest possible start point */
   int                   duration,           /**< duration of the job */
   int                   demand,             /**< demand of the job */
   SCIP_Bool*            infeasible          /**< pointer store if the job cannot be scheduled */
   )
{
   int starttime;
   int pos;

   assert(profile != NULL);
   assert(lb >= 0);
   assert(lb <= ub);
   assert(duration >= 0);
   assert(demand >= 0);
   assert(infeasible != NULL);
   assert(profile->timepoints[profile->ntimepoints-1] > ub);

   if( duration == 0 || demand == 0 )
      return ub;

   starttime = ub;
   (void)SCIPprofileFindLowerBound(profile, starttime, &pos);
   assert(profile->timepoints[pos] <= starttime);

   (*infeasible) = TRUE;

   while( (*infeasible) && starttime >= lb )
   {
      if( SCIPprofileIsFeasibleStart(profile, starttime, duration, demand, &pos) )
      {
         (*infeasible) = FALSE;
         return starttime;
      }
      assert(pos >= 0);

      /* the job did not fit into the profile since at time point "pos" not enough capacity is available;
       * therefore we can proceed with the next time point  */
      assert(profile->freecapacities[pos] < demand);

      starttime = profile->timepoints[pos] - duration;
   }

   assert(*infeasible || starttime >= lb);

   return starttime;
}

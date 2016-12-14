/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   rbtree.c
 * @brief  intrusive red black tree datastructure
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/rbtree.h"

#define RED              ((uintptr_t)0x1u)
#define BLACK            ((uintptr_t)0x0u)
#define COLOR(node)      ((node)->parent & RED)
#define IS_RED(node)     ( (node) != NULL && COLOR(node) )
#define IS_BLACK(node)   ( (node) == NULL || !COLOR(node) )
#define MAKE_RED(node)   do { (node)->parent |= RED; } while(0)
#define MAKE_BLACK(node) do { (node)->parent &= ~RED; } while(0)
#define LEFT             0
#define RIGHT            1
#define OPPOSITE(dir)    ( 1 - (dir) )
#define PARENT(node)     ( (SCIP_RBTREENODE*)((node)->parent & ~RED) )
#define SET_PARENT(n, p) do { (n)->parent = (uintptr_t)(p) | COLOR(n); } while(0)
#define SET_COLOR(n, c)  do { if( c == RED ) { MAKE_RED(n); } else { MAKE_BLACK(n); } } while(0)


static
void rbRotate(
   SCIP_RBTREENODE**     root,
   SCIP_RBTREENODE*      x,
   int                   dir
   )
{
   SCIP_RBTREENODE* p;
   SCIP_RBTREENODE* y = x->child[OPPOSITE(dir)];
   x->child[OPPOSITE(dir)] = y->child[dir];
   if( y->child[dir] != NULL )
   {
      SET_PARENT(y->child[dir], x);
   }

   p = PARENT(x);
   SET_PARENT(y, p);

   if( p == NULL )
      *root = y;
   else if( x == p->child[dir] )
      p->child[dir] = y;
   else
      p->child[OPPOSITE(dir)] = y;

   y->child[dir] = x;
   SET_PARENT(x, y);
}

static
void rbInsertFixup(
   SCIP_RBTREENODE**     root,
   SCIP_RBTREENODE*      z
   )
{
   SCIP_RBTREENODE* p;
   p = PARENT(z);

   while( IS_RED(p) )
   {
      SCIP_RBTREENODE* pp;
      SCIP_RBTREENODE* y;
      int dir;

      pp = PARENT(p);
      dir = p == pp->child[LEFT] ? RIGHT : LEFT;

      y = pp->child[dir];
      if( IS_RED(y) )
      {
         MAKE_BLACK(p);
         MAKE_BLACK(y);
         MAKE_RED(pp);
         z = pp;
      }
      else
      {
         if( z == p->child[dir] )
         {
            z = p;
            rbRotate(root, z, OPPOSITE(dir));
            p = PARENT(z);
            pp = PARENT(p);
         }

         MAKE_BLACK(p);
         MAKE_RED(pp);
         rbRotate(root, pp, dir);
      }

      p = PARENT(z);
   }

   MAKE_BLACK(*root);
}

static
void rbDeleteFixup(
   SCIP_RBTREENODE**     root,
   SCIP_RBTREENODE*      x,
   SCIP_RBTREENODE*      nil
   )
{
   while( x != *root && IS_BLACK(x) )
   {
      SCIP_RBTREENODE* p;
      SCIP_RBTREENODE* w;
      int dir;

      p = PARENT(x == NULL ? nil : x);
      dir = x == p->child[LEFT] ? RIGHT : LEFT;

      w = p->child[dir];
      if( IS_RED(w) )
      {
         MAKE_BLACK(w);
         MAKE_RED(p);
         rbRotate(root, p, OPPOSITE(dir));
         assert(p == PARENT(x == NULL ? nil : x));
         w = p->child[dir];
      }

      if( IS_BLACK(w->child[LEFT]) && IS_BLACK(w->child[RIGHT]) )
      {
         MAKE_RED(w);
         x = p;
      }
      else
      {
         if( IS_BLACK(w->child[dir]) )
         {
            MAKE_BLACK(w->child[OPPOSITE(dir)]);
            MAKE_RED(w);
            rbRotate(root, w, dir);
            assert(p == PARENT(x == NULL ? nil : x));
            w = p->child[dir];
         }
         SET_COLOR(w, COLOR(p));
         MAKE_BLACK(p);
         MAKE_BLACK(w->child[dir]);
         rbRotate(root, p, OPPOSITE(dir));
         x = *root;
      }
   }

   if( x != NULL )
   {
      MAKE_BLACK(x);
   }
}

static
void rbTransplant(
   SCIP_RBTREENODE**     root,             /**< chunk block */
   SCIP_RBTREENODE*      u,
   SCIP_RBTREENODE*      v,
   SCIP_RBTREENODE*      nil
   )
{
   SCIP_RBTREENODE* up;

   up = PARENT(u);

   if( up == NULL )
      *root = v;
   else if( u == up->child[LEFT] )
      up->child[LEFT] = v;
   else
      up->child[RIGHT] = v;

   if( v == NULL )
      v = nil;

   SET_PARENT(v, up);
}

SCIP_RBTREENODE* SCIPrbtreeFirst_call(
   SCIP_RBTREENODE*      root
   )
{
   if( root == NULL )
      return NULL;

   while(root->child[LEFT] != NULL)
      root = root->child[LEFT];

   return root;
}

SCIP_RBTREENODE* SCIPrbtreeLast_call(
   SCIP_RBTREENODE*      root
   )
{
   if( root == NULL )
      return NULL;

   while(root->child[RIGHT] != NULL)
      root = root->child[RIGHT];

   return root;
}

SCIP_RBTREENODE* SCIPrbtreeSuccessor_call(
   SCIP_RBTREENODE* x
   )
{
   SCIP_RBTREENODE* y;
   if( x->child[RIGHT] != NULL )
      return SCIPrbtreeFirst_call(x->child[RIGHT]);

   y = PARENT(x);

   while( y != NULL && x == y->child[RIGHT] )
   {
      x = y;
      y = PARENT(y);
   }

   return y;
}

SCIP_RBTREENODE* SCIPrbtreePredecessor_call(
   SCIP_RBTREENODE* x
   )
{
   SCIP_RBTREENODE* y;
   if( x->child[LEFT] != NULL )
      return SCIPrbtreeLast_call(x->child[LEFT]);

   y = PARENT(x);

   while( y != NULL && x == y->child[LEFT] )
   {
      x = y;
      y = PARENT(y);
   }

   return y;
}

void SCIPrbtreeDelete_call(
   SCIP_RBTREENODE**     root,
   SCIP_RBTREENODE*      node
   )
{
   SCIP_RBTREENODE nil;
   SCIP_RBTREENODE* y;
   SCIP_RBTREENODE* x;
   unsigned int yorigcolor;

   y = node;
   yorigcolor = COLOR(y);

   if( node->child[LEFT] == NULL )
   {
      x = node->child[RIGHT];
      rbTransplant(root, node, x, &nil);
   }
   else if( node->child[RIGHT] == NULL )
   {
      x = node->child[LEFT];
      rbTransplant(root, node, x, &nil);
   }
   else
   {
      y = SCIPrbtreeFirst(node->child[RIGHT]);
      yorigcolor = COLOR(y);
      x = y->child[RIGHT];
      if( PARENT(y) == node )
      {
         SET_PARENT(x == NULL ? &nil : x, y);
      }
      else
      {
         rbTransplant(root, y, y->child[RIGHT], &nil);
         y->child[RIGHT] = node->child[RIGHT];
         SET_PARENT(y->child[RIGHT], y);
      }
      rbTransplant(root, node, y, &nil);
      y->child[LEFT] = node->child[LEFT];
      SET_PARENT(y->child[LEFT], y);
      SET_COLOR(y, COLOR(node));
   }

   if( yorigcolor == BLACK )
      rbDeleteFixup(root, x, &nil);
}

void SCIPrbtreeInsert_call(
   SCIP_RBTREENODE**     root,
   SCIP_RBTREENODE*      parent,
   int                   pos,
   SCIP_RBTREENODE*      node
   )
{
   SET_PARENT(node, parent);
   if( parent == NULL )
      *root = node;
   else if( pos > 0 )
      parent->child[LEFT] = node;
   else
      parent->child[RIGHT] = node;

   node->child[LEFT] = NULL;
   node->child[RIGHT] = NULL;
   MAKE_RED(node);
   rbInsertFixup(root, node);
}

#define SORTTPL_NAMEEXT Int_call
#define SORTTPL_KEYTYPE int
#include "rbtreetpl.c"

#define SORTTPL_NAMEEXT Real_call
#define SORTTPL_KEYTYPE SCIP_Real
#include "rbtreetpl.c"

#define SORTTPL_NAMEEXT Ptr_call
#define SORTTPL_KEYTYPE void*
#define SORTTPL_PTRCOMP
#include "rbtreetpl.c"

#define SORTTPL_NAMEEXT Elem_call
#define RBTREE_NO_KEY
#define SORTTPL_PTRCOMP
#include "rbtreetpl.c"

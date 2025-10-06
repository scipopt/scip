/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file pqueue.h
 * @brief  class for priority queues
 * @author Andreas Bley
 * @author Marc Pfetsch
 */

#ifndef _PQUEUE_H
#define _PQUEUE_H

#include <algorithm>
#include <functional>

namespace std
{
   ///
   template<typename Key,
            typename Data,
            typename Compare = less<Key> >

   class pqueue /*lint --e{3713}*/
   {
      private:

      //--------------------
      // item node class
      //--------------------
      class node /*lint --e{3713}*/
      {
         friend class pqueue;

      public:
         //
         node(
            const Key& k,
            const Data& d ):
            key   (k),
            data  (d),
            sleft (0),
            sright(0),
            left  (NULL),
            right (NULL),
            father(NULL)
         {}

         //
         ~node()
         {
            left = nullptr;
            right = nullptr;
            father = nullptr;
         }

      private:
         //
         void delete_children_recursive()
         {
            if ( left != NULL )
            {
               left->delete_children_recursive();
               delete left;
               left = NULL;
            }
            if ( right != NULL )
            {
               right->delete_children_recursive();
               delete right;
               right = NULL;
            }
         }

         Key   key;
         Data  data;
         int   sleft;
         int   sright;
         node* left;
         node* right;
         node* father;
      };

   public:

      typedef node* pqueue_item;

   private:

      node*   root;
      Compare compare;

   public:

      /** Default constructor, creates empty priority queue. */
      pqueue():
         root( NULL )
      {} /*lint !e1401*/

      /** Destructs queue */
      ~pqueue()
      {
         clear(); /*lint !e1551*/
      } /*lint !e1579*/

      /** Empties queue */
      void clear()
      {
         if ( root != NULL )
         {
            root->delete_children_recursive();
            delete root;
            root = NULL;
         }
      }

      /** Returns true if the pqueue is empty. */
      bool empty() const
      {
         return ( root == NULL );
      }

      /** Returns size of queue. */
      int size() const
      {
         return ( root == NULL ? 0 : root->sleft + root->sright + 1 );
      }

      /** Returns key of queue item. */
      const Key& get_key(
         pqueue_item it
         ) const
      {
         assert( it != NULL );
         return it->key;
      }

      /** Returns data of queue item. */
      const Data& get_data(
         pqueue_item it
         ) const
      {
         assert( it != NULL );
         return it->data;
      }

      /** Returns queue item at top (with lowers key). */
      pqueue_item top()
      {
         return root; /*lint !e1535*//*lint !e1806*/
      }

      /** Inserts a new entry into the queue, returns new item */
      pqueue_item  insert(
         const Key&  key,
         const Data& data
         )
      {
         node* nn = NULL;
         if ( root == NULL )
         {
            nn = new node(key,data);
            if ( nn == NULL ) /*lint !e774*/
               throw std::bad_alloc();
            root = nn;
         }
         else
            nn = create_new_node(key, data, root);

         rotate_backward(nn);
         return nn;
      }

      /** Reduces the key a queue item. */
      void decrease_key(
         pqueue_item item,
         const Key&  new_key
         )
      {
         assert( item );
         assert( compare(new_key, item->key) );

         item->key = new_key;
         rotate_backward(item);
      }

      /** Removes the topmost item from the queue. */
      void pop()
      {
         assert ( root != NULL );
         remove( root );
      }

      /** Removes the item from the queue */
      void remove(
         node* item
         )
      {
         assert ( item != NULL );
         assert ( root != NULL );

         bool goto_left  = ( item->left  != NULL );
         bool goto_right = ( item->right != NULL );
         if ( goto_left && goto_right )
         {
            goto_right = ( compare( item->right->key, item->left->key ) );
            goto_left  = ! goto_right;
         }
         if ( goto_right )
         {
            swap_with_father( item->right );
            remove( item );
            return;
         }
         if ( goto_left )
         {
            swap_with_father( item->left );
            remove( item );
            return;
         }
         // at leave: remove and update all sizes
         for (node* n = item, *f = n->father; f != NULL; n = f, f = n->father) /*lint !e440*/ /*lint !e2840*/
         {
            if ( f->left == n )
            {
               f->sleft -= 1;
            }
            else
            {
               assert( f->right == n );
               f->sright -= 1;
            }
         }
         if ( item->father )
         {
            if ( item->father->left == item )
            {
               assert( item->father->sleft == 0 );
               item->father->left = NULL;
            }
            else
            {
               assert( item->father->right == item );
               assert( item->father->sright == 0 );
               item->father->right = NULL;
            }
         }
         else
         {
            assert( item == root );
            root = NULL;
         }
         delete item;
      }


   private:

      /** creates new element in the tree such that tree remains balanced */
      node* create_new_node(
         const Key&  key,
         const Data& data,
         node*       subproblem
         )
      {
         assert( subproblem != NULL );

         if ( subproblem->sleft == 0 )
         {
            assert( subproblem->left == NULL );

            node* nn = new node(key,data);
            subproblem->left  = nn;
            subproblem->sleft = 1;
            nn->father        = subproblem;
            return nn;
         }
         if ( subproblem->sright == 0 )
         {
            assert( subproblem->right == NULL );

            node* nn = new node(key,data);
            subproblem->right  = nn;
            subproblem->sright = 1;
            nn->father         = subproblem;
            return nn;
         }
         assert( subproblem->left  != NULL );
         assert( subproblem->right != NULL );

         if ( subproblem->sleft <= subproblem->sright )
         {
            subproblem->sleft += 1;
            return create_new_node(key, data, subproblem->left);
         }

         subproblem->sright += 1;
         return create_new_node(key, data, subproblem->right);
      }


      void swap_with_father(
         node* n1
         )
      {
         int   n1_sleft  = n1->sleft;
         int   n1_sright = n1->sright;
         node* n1_left   = n1->left;
         node* n1_right  = n1->right;
         node* n1_father = n1->father;
         assert( n1_father != NULL );
         assert( n1_father->left == n1 || n1_father->right == n1 );

         if ( root == n1_father )
            root = n1;

         if ( n1_father->left == n1 )
         {
            n1->left   = n1_father;
            n1->right  = n1_father->right;
         }
         else
         {
            assert( n1_father->right == n1 );

            n1->left          = n1_father->left;
            n1->right         = n1_father;
         }
         n1_father->left   = n1_left;
         n1_father->right  = n1_right;

         n1->sleft         = n1_father->sleft;
         n1->sright        = n1_father->sright;
         n1_father->sleft  = n1_sleft;
         n1_father->sright = n1_sright;

         n1->father        = n1_father->father;
         n1_father->father = n1;

         if ( n1->left )
            n1->left-> father = n1;
         if ( n1->right )
            n1->right->father = n1;
         if ( n1_father->left  )
            n1_father->left->father  = n1_father;
         if ( n1_father->right )
            n1_father->right->father = n1_father;
         if ( n1->father )
         {
            if ( n1->father->left == n1_father )
               n1->father->left = n1;
            if ( n1->father->right == n1_father )
               n1->father->right = n1;
         }
      }

      void rotate_backward(
         node* item
         )
      {
         assert( item != NULL );

         if ( item->father )
         {
            if ( ! compare( item->father->key, item->key ) )
            {
               swap_with_father( item );
               rotate_backward( item );
            }
         }
      }
   }; /*lint !e1712*/

} // namespace std

#endif /* _PQUEUE_H */

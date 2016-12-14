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

#ifndef SORTTPL_NAMEEXT
#error You need to define SORTTPL_NAMEEXT.
#endif
#if !defined(SORTTPL_KEYTYPE) && !defined(RBTREE_NO_KEY)
#error You need to define SORTTPL_KEYTYPE or RBTREE_NO_KEY.
#endif

#ifdef SORTTPL_EXPANDNAME
#undef SORTTPL_EXPANDNAME
#endif
#ifdef SORTTPL_NAME
#undef SORTTPL_NAME
#endif

/* enabling and disabling additional lines in the code */
#ifdef SORTTPL_PTRCOMP
#define SORTTPL_HASPTRCOMP(x)    x
#define SORTTPL_HASPTRCOMPPAR(x) x,
#else
#define SORTTPL_HASPTRCOMP(x)    /**/
#define SORTTPL_HASPTRCOMPPAR(x) /**/
#endif
#ifdef SORTTPL_INDCOMP
#define SORTTPL_HASINDCOMP(x)    x
#define SORTTPL_HASINDCOMPPAR(x) x,
#else
#define SORTTPL_HASINDCOMP(x)    /**/
#define SORTTPL_HASINDCOMPPAR(x) /**/
#endif


/* the two-step macro definition is needed, such that macro arguments
 * get expanded by prescan of the C preprocessor (see "info cpp",
 * chapter 3.10.6: Argument Prescan)
 */
#define SORTTPL_EXPANDNAME(method, methodname) \
   method ## methodname
#define SORTTPL_NAME(method, methodname) \
  SORTTPL_EXPANDNAME(method, methodname)

/* comparator method */
#ifdef SORTTPL_PTRCOMP
#ifdef SORTTPL_BACKWARDS
#define SORTTPL_ISBETTER(x,y) (ptrcomp((x), (y)) > 0)
#define SORTTPL_ISWORSE(x,y) (ptrcomp((x), (y)) < 0)
#else
#define SORTTPL_ISBETTER(x,y) (ptrcomp((x), (y)) < 0)
#define SORTTPL_ISWORSE(x,y) (ptrcomp((x), (y)) > 0)
#endif
#else
#ifdef SORTTPL_INDCOMP
#ifdef SORTTPL_BACKWARDS
#define SORTTPL_ISBETTER(x,y) (indcomp(dataptr, (x), (y)) > 0)
#define SORTTPL_ISWORSE(x,y) (indcomp(dataptr, (x), (y)) < 0)
#else
#define SORTTPL_ISBETTER(x,y) (indcomp(dataptr, (x), (y)) < 0)
#define SORTTPL_ISWORSE(x,y) (indcomp(dataptr, (x), (y)) > 0)
#endif
#else
#ifdef SORTTPL_BACKWARDS
#define SORTTPL_ISBETTER(x,y) ((x) > (y))
#define SORTTPL_ISWORSE(x,y) ((x) < (y))
#else
#define SORTTPL_ISBETTER(x,y) ((x) < (y))
#define SORTTPL_ISWORSE(x,y) ((x) > (y))
#endif
#endif
#endif


#ifdef RBTREENODE_KEY
#undef RBTREENODE_KEY
#endif

#ifdef RBTREE_NO_KEY

#define RBTREENODE_KEY(n) ((void*)(n))
#define SORTTPL_KEYTYPE SCIP_RBTREENODE*

#else

struct SORTTPL_NAME(__SCIP_RBTREENODE_, SORTTPL_NAMEEXT)
{
   SCIP_RBTREE_KEY(SORTTPL_KEYTYPE key);
};

#define RBTREENODE_KEY(n) ( ((struct SORTTPL_NAME(__SCIP_RBTREENODE_, SORTTPL_NAMEEXT)*)(n))->key )
#endif


int SORTTPL_NAME(SCIPrbtreeFind, SORTTPL_NAMEEXT)(
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   SCIP_RBTREENODE*      root,
   SORTTPL_KEYTYPE       key,
   SCIP_RBTREENODE**     node
   )
{
   SCIP_RBTREENODE* x;

   *node = NULL;
   x = root;

   while( x != NULL )
   {
      *node = x;
      if( SORTTPL_ISBETTER(key, RBTREENODE_KEY(x)) )
         x = x->child[LEFT];
      else if( SORTTPL_ISWORSE(key, RBTREENODE_KEY(x)) )
         x = x->child[RIGHT];
      else
         return 0;
   }

   if( *node != NULL && SORTTPL_ISBETTER(key, RBTREENODE_KEY(*node)) )
      return 1;

   return -1;
}

/* undefine template parameters and local defines */
#undef SORTTPL_NAMEEXT
#undef SORTTPL_KEYTYPE
#undef SORTTPL_PTRCOMP
#undef SORTTPL_INDCOMP
#undef SORTTPL_HASPTRCOMP
#undef SORTTPL_HASINDCOMP
#undef SORTTPL_HASPTRCOMPPAR
#undef SORTTPL_HASINDCOMPPAR
#undef SORTTPL_ISBETTER
#undef SORTTPL_ISWORSE
#undef SORTTPL_BACKWARDS
#undef RBTREENODE_KEY

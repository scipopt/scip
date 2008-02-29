/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sorttpl.c,v 1.1 2008/02/29 11:03:50 bzfpfend Exp $"

/**@file   sorttpl.c
 * @brief  template functions for sorting
 * @author Michael Winkler
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* template parameters that have to be passed in as #define's:
 * #define SORTTPL_METHOD       <name>     name of the SCIP method that should be generated
 * #define SORTTPL_KEYTYPE      <type>     data type of the key array
 * #define SORTTPL_FIELD1TYPE   <type>     data type of first additional array which should be sorted in the same way (optional)
 * #define SORTTPL_FIELD2TYPE   <type>     data type of second additional array which should be sorted in the same way (optional)
 * #define SORTTPL_FIELD3TYPE   <type>     data type of third additional array which should be sorted in the same way (optional)
 * #define SORTTPL_PTRCOMP                 ptrcomp method should be used for comparisons (optional)
 * #define SORTTPL_INDCOMP                 indcomp method should be used for comparisons (optional)
 */

#define SORTTPL_SHELLSORTMAX 25

#ifndef SORTTPL_METHOD
#error You need to define the SORTTPL_METHOD name.
#endif
#ifndef SORTTPL_KEYTYPE
#error You need to define the SORTTPL_KEYTYPE.
#endif

#ifdef SORTTPL_EXPANDNAME
#undef SORTTPL_EXPANDNAME
#endif
#ifdef SORTTPL_NAME
#undef SORTTPL_NAME
#endif

/* enabling and disabling additional lines in the code */
#ifdef SORTTPL_FIELD1TYPE
#define SORTTPL_HASFIELD1(x)    x
#define SORTTPL_HASFIELD1PAR(x) x,
#else
#define SORTTPL_HASFIELD1(x)    /**/
#define SORTTPL_HASFIELD1PAR(x) /**/
#define SORTTPL_FIELD1NAME undef
#endif
#ifdef SORTTPL_FIELD2TYPE
#define SORTTPL_HASFIELD2(x)    x
#define SORTTPL_HASFIELD2PAR(x) x,
#else
#define SORTTPL_HASFIELD2(x)    /**/
#define SORTTPL_HASFIELD2PAR(x) /**/
#define SORTTPL_FIELD2NAME undef
#endif
#ifdef SORTTPL_FIELD3TYPE
#define SORTTPL_HASFIELD3(x)    x
#define SORTTPL_HASFIELD3PAR(x) x,
#else
#define SORTTPL_HASFIELD3(x)    /**/
#define SORTTPL_HASFIELD3PAR(x) /**/
#define SORTTPL_FIELD3NAME undef
#endif
#ifdef SORTTPL_PTRCOMP
#define SORTTPL_HASPTRCOMP(x)    x
#define SORTTPL_HASPTRCOMPPAR(x) x,
#define SORTTPL_PTRCOMPNAME Yes
#else
#define SORTTPL_HASPTRCOMP(x)    /**/
#define SORTTPL_HASPTRCOMPPAR(x) /**/
#define SORTTPL_PTRCOMPNAME No
#endif
#ifdef SORTTPL_INDCOMP
#define SORTTPL_HASINDCOMP(x)    x
#define SORTTPL_HASINDCOMPPAR(x) x,
#define SORTTPL_INDCOMPNAME Yes
#else
#define SORTTPL_HASINDCOMP(x)    /**/
#define SORTTPL_HASINDCOMPPAR(x) /**/
#define SORTTPL_INDCOMPNAME No
#endif

/* names of the fields for method names */
#ifndef SORTTPL_KEYNAME
#define SORTTPL_KEYNAME _ ## SORTTPL_KEYTYPE
#endif
#ifndef SORTTPL_FIELD1NAME
#define SORTTPL_FIELD1NAME _ ## SORTTPL_FIELD1TYPE
#endif
#ifndef SORTTPL_FIELD2NAME
#define SORTTPL_FIELD2NAME _ ## SORTTPL_FIELD2TYPE
#endif
#ifndef SORTTPL_FIELD3NAME
#define SORTTPL_FIELD3NAME _ ## SORTTPL_FIELD3TYPE
#endif


/* the two-step macro definition is needed, such that macro arguments
 * get expanded by prescan of the C preprocessor (see "info cpp",
 * chapter 3.10.6: Argument Prescan)
 */
#define SORTTPL_EXPANDNAME(method, keyname, field1name, field2name, field3name, ptrcomp, indcomp) \
   sorttpl_ ## method ## _ ## keyname ## _ ## field1name ## _ ## field2name ## _ ## field3name ## _ ## ptrcomp
#define SORTTPL_NAME(method, keyname, field1name, field2name, field3name, ptrcomp, indcomp) \
   SORTTPL_EXPANDNAME(method, keyname, field1name, field2name, field3name, ptrcomp, indcomp)

/* comparator method */
#ifdef SORTTPL_PTRCOMP
#define SORTTPL_ISSMALLER(x,y) (ptrcomp((x), (y)) < 0)
#else
#ifdef SORTTPL_INDCOMP
#define SORTTPL_ISSMALLER(x,y) (indcomp(dataptr, (x), (y)) < 0)
#else
#define SORTTPL_ISSMALLER(x,y) ((x) < (y))
#endif
#endif

/* swapping two variables */
#define SORTTPL_SWAP(T,x,y) \
   {                \
      T temp = x;   \
      x = y;        \
      y = temp;     \
   }



/** shellsort an array of data elements; use it only for arrays smaller than 25 entries */
static
void SORTTPL_NAME(shellSort, SORTTPL_KEYNAME, SORTTPL_FIELD1NAME, SORTTPL_FIELD2NAME, SORTTPL_FIELD3NAME, SORTTPL_PTRCOMPNAME, SORTTPL_INDCOMPNAME)
(
   SORTTPL_KEYTYPE*      key,                /**< pointer to data array that defines the order */
   SORTTPL_HASFIELD1PAR(  SORTTPL_FIELD1TYPE*    field1 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD2PAR(  SORTTPL_FIELD2TYPE*    field2 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD3PAR(  SORTTPL_FIELD3TYPE*    field3 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   int                   start,              /**< starting index */
   int                   end                 /**< ending index */
   )
{
   static const int incs[3] = {1, 5, 19}; /* sequence of increments */
   int k;

   assert(start <= end);

   for( k = 3; k >= 0; --k )
   {
      int h;
      int i;

      for( h = incs[k], i = h + start; i <= end; ++i )
      {
         SORTTPL_KEYTYPE tempkey = key[i];
         SORTTPL_HASFIELD1( SORTTPL_FIELD1TYPE tempfield1 = field1[i]; );
         SORTTPL_HASFIELD2( SORTTPL_FIELD2TYPE tempfield2 = field2[i]; );
         SORTTPL_HASFIELD3( SORTTPL_FIELD3TYPE tempfield3 = field3[i]; );
         int j;

         j = i;
         while( j >= h && SORTTPL_ISSMALLER(tempkey, key[j-h]) )
         {
            key[j] = key[j-h];
            SORTTPL_HASFIELD1( field1[j] = field1[j-h]; );
            SORTTPL_HASFIELD2( field2[j] = field2[j-h]; );
            SORTTPL_HASFIELD3( field3[j] = field3[j-h]; );
            j -= h;
         }
          
         key[j] = tempkey;
         SORTTPL_HASFIELD1( field1[j] = tempfield1; );
         SORTTPL_HASFIELD2( field2[j] = tempfield2; );
         SORTTPL_HASFIELD3( field3[j] = tempfield3; );
      }
   }
}


/** quicksort an array of pointers; pivot is the medial element */
static
void SORTTPL_NAME(qSort, SORTTPL_KEYNAME, SORTTPL_FIELD1NAME, SORTTPL_FIELD2NAME, SORTTPL_FIELD3NAME, SORTTPL_PTRCOMPNAME, SORTTPL_INDCOMPNAME)
(
   SORTTPL_KEYTYPE*      key,                /**< pointer to data array that defines the order */
   SORTTPL_HASFIELD1PAR(  SORTTPL_FIELD1TYPE*    field1 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD2PAR(  SORTTPL_FIELD2TYPE*    field2 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD3PAR(  SORTTPL_FIELD3TYPE*    field3 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   int                   start,              /**< starting index */
   int                   end                 /**< ending index */
   )
{
   assert(start <= end);

   /* use quick sort for long lists */
   while( end - start >= SORTTPL_SHELLSORTMAX )
   {
      SORTTPL_KEYTYPE pivotkey;
      int lo;
      int hi;
      int mid;

      /* select pivot element */
      mid = (start+end)/2;
      pivotkey = key[mid];

      /* partition the array into elements < pivot [start,hi] and elements >= pivot [lo,end] */
      lo = start;
      hi = end;
      for( ;; )
      {
         while( lo <= end && SORTTPL_ISSMALLER(key[lo], pivotkey) )
            lo++;
         while( hi >= start && !SORTTPL_ISSMALLER(key[hi], pivotkey) )
            hi--;
      
         if( lo >= hi )
            break;

         SORTTPL_SWAP(SORTTPL_KEYTYPE, key[lo], key[hi]);
         SORTTPL_HASFIELD1( SORTTPL_SWAP(SORTTPL_FIELD1TYPE, field1[lo], field1[hi]) );
         SORTTPL_HASFIELD2( SORTTPL_SWAP(SORTTPL_FIELD2TYPE, field2[lo], field2[hi]) );
         SORTTPL_HASFIELD3( SORTTPL_SWAP(SORTTPL_FIELD3TYPE, field3[lo], field3[hi]) );

         lo++;
         hi--;
      }
      assert(hi == lo-1);

      /* skip entries which are equal to the pivot element (three partitions, <, =, > than pivot)*/
      while( lo <= end && !SORTTPL_ISSMALLER(pivotkey, key[lo]) )
         lo++;

      /* make sure that we have at least one element in the smaller partition */
      if( lo == start )
      {
         /* everything is greater or equal than the pivot element: move pivot to the left (degenerate case) */
         assert(!SORTTPL_ISSMALLER(key[mid], pivotkey)); /* the pivot element did not change its position */
         assert(!SORTTPL_ISSMALLER(pivotkey, key[mid]));
         SORTTPL_SWAP(SORTTPL_KEYTYPE, key[lo], key[mid]);
         SORTTPL_HASFIELD1( SORTTPL_SWAP(SORTTPL_FIELD1TYPE, field1[lo], field1[mid]) );
         SORTTPL_HASFIELD2( SORTTPL_SWAP(SORTTPL_FIELD2TYPE, field2[lo], field2[mid]) );
         SORTTPL_HASFIELD3( SORTTPL_SWAP(SORTTPL_FIELD3TYPE, field3[lo], field3[mid]) );
         lo++;
      }

      /* sort the smaller partition by a recursive call, sort the larger part without recursion */
      if( hi - start <= end - lo )
      {
         /* sort [start,hi] with a recursive call */
         if( start < hi )
         {
            SORTTPL_NAME(qSort, SORTTPL_KEYNAME, SORTTPL_FIELD1NAME, SORTTPL_FIELD2NAME, SORTTPL_FIELD3NAME, SORTTPL_PTRCOMPNAME, SORTTPL_INDCOMPNAME)
               (key,
                SORTTPL_HASFIELD1PAR(field1)
                SORTTPL_HASFIELD2PAR(field2)
                SORTTPL_HASFIELD3PAR(field3)
                SORTTPL_HASPTRCOMPPAR(ptrcomp)
                SORTTPL_HASINDCOMPPAR(indcomp)
                SORTTPL_HASINDCOMPPAR(dataptr)
                start, hi);
         }

         /* now focus on the larger part [lo,end] */
         start = lo;
      }
      else
      {
         if( lo < end )
         {
            /* sort [lo,end] with a recursive call */
            SORTTPL_NAME(qSort, SORTTPL_KEYNAME, SORTTPL_FIELD1NAME, SORTTPL_FIELD2NAME, SORTTPL_FIELD3NAME, SORTTPL_PTRCOMPNAME, SORTTPL_INDCOMPNAME)
               (key,
                SORTTPL_HASFIELD1PAR(field1)
                SORTTPL_HASFIELD2PAR(field2)
                SORTTPL_HASFIELD3PAR(field3)
                SORTTPL_HASPTRCOMPPAR(ptrcomp)
                SORTTPL_HASINDCOMPPAR(indcomp)
                SORTTPL_HASINDCOMPPAR(dataptr)
                lo, end);
         }

         /* now focus on the larger part [start,hi] */
         end = hi;
      }
   }

   /* use shell sort on the remaining small list */
   if( end - start >= 1 )
   {
      SORTTPL_NAME(shellSort, SORTTPL_KEYNAME, SORTTPL_FIELD1NAME, SORTTPL_FIELD2NAME, SORTTPL_FIELD3NAME, SORTTPL_PTRCOMPNAME, SORTTPL_INDCOMPNAME)
         (key,
          SORTTPL_HASFIELD1PAR(field1)
          SORTTPL_HASFIELD2PAR(field2)
          SORTTPL_HASFIELD3PAR(field3)
          SORTTPL_HASPTRCOMPPAR(ptrcomp)
          SORTTPL_HASINDCOMPPAR(indcomp)
          SORTTPL_HASINDCOMPPAR(dataptr)
          start, end);
   }
}


/** verifies that an array is indeed sorted */
static
void SORTTPL_NAME(checkSort, SORTTPL_KEYNAME, SORTTPL_FIELD1NAME, SORTTPL_FIELD2NAME, SORTTPL_FIELD3NAME, SORTTPL_PTRCOMPNAME, SORTTPL_INDCOMPNAME)
(
   SORTTPL_KEYTYPE*      key,                /**< pointer to data array that defines the order */
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< length of the array */
   )
{
#ifndef NDEBUG
   int i;

   for( i = 0; i < len-1; i++ )
   {
      assert(!SORTTPL_ISSMALLER(key[i+1], key[i]));
   }
#endif
}

/** sorts array 'key' and performs the same permutations on the additional 'field' arrays */
void SORTTPL_METHOD (
   SORTTPL_KEYTYPE*      key,                /**< pointer to data array that defines the order */
   SORTTPL_HASFIELD1PAR(  SORTTPL_FIELD1TYPE*    field1 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD2PAR(  SORTTPL_FIELD2TYPE*    field2 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD3PAR(  SORTTPL_FIELD3TYPE*    field3 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< length of arrays */
   )
{
   SORTTPL_NAME(qSort, SORTTPL_KEYNAME, SORTTPL_FIELD1NAME, SORTTPL_FIELD2NAME, SORTTPL_FIELD3NAME, SORTTPL_PTRCOMPNAME, SORTTPL_INDCOMPNAME)
      (key,
       SORTTPL_HASFIELD1PAR(field1)
       SORTTPL_HASFIELD2PAR(field2)
       SORTTPL_HASFIELD3PAR(field3)
       SORTTPL_HASPTRCOMPPAR(ptrcomp)
       SORTTPL_HASINDCOMPPAR(indcomp)
       SORTTPL_HASINDCOMPPAR(dataptr)
       0, len-1);
   SORTTPL_NAME(checkSort, SORTTPL_KEYNAME, SORTTPL_FIELD1NAME, SORTTPL_FIELD2NAME, SORTTPL_FIELD3NAME, SORTTPL_PTRCOMPNAME, SORTTPL_INDCOMPNAME)
      (key,
       SORTTPL_HASPTRCOMPPAR(ptrcomp)
       SORTTPL_HASINDCOMPPAR(indcomp)
       SORTTPL_HASINDCOMPPAR(dataptr)
       len);
}


/* undefine template parameters and local defines */
#undef SORTTPL_METHOD
#undef SORTTPL_KEYTYPE
#undef SORTTPL_FIELD1TYPE
#undef SORTTPL_FIELD2TYPE
#undef SORTTPL_FIELD3TYPE
#undef SORTTPL_PTRCOMP
#undef SORTTPL_INDCOMP
#undef SORTTPL_KEYNAME
#undef SORTTPL_FIELD1NAME
#undef SORTTPL_FIELD2NAME
#undef SORTTPL_FIELD3NAME
#undef SORTTPL_PTRCOMPNAME
#undef SORTTPL_INDCOMPNAME
#undef SORTTPL_HASFIELD1
#undef SORTTPL_HASFIELD2
#undef SORTTPL_HASFIELD3
#undef SORTTPL_HASPTRCOMP
#undef SORTTPL_HASINDCOMP
#undef SORTTPL_HASFIELD1PAR
#undef SORTTPL_HASFIELD2PAR
#undef SORTTPL_HASFIELD3PAR
#undef SORTTPL_HASPTRCOMPPAR
#undef SORTTPL_HASINDCOMPPAR
#undef SORTTPL_ISSMALLER
#undef SORTTPL_SWAP
#undef SORTTPL_SHELLSORTMAX

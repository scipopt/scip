/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* prints short statistics (callback, preprocessing, adding cuts) */
/* // #define SCIP_DEBUG */
/* // #define ZEROHALF__PRINT_STATISTICS */ /**< print statistics */

/**
 * @file   sepa_zerohalf.c
 * @brief  {0,1/2}-cuts separator
 * @author Manuel Kutschka
 * @author Kati Wolter
 *
 * {0,1/2}-Chv'atal-Gomory cuts separator. It solves the following separation problem:
 *
 *
 * Given an integer program 
 *
 * - min { c^T x : Ax <= b, x >= 0, x integer }
 *
 * and a fractional solution x* of its LP relaxation.
 *
 * Find a weightvector u whose entries u_i are either 0 or 1/2 such that the following inequality is valid for all integral solutions and violated by x*
 *
 * - floor(u^T A) x <= floor(u^T b)
 *
 * or (if exact methods are used) give a proof that no such inequality exists
 *
 *
 *
 * References:
 * - Alberto Caprara, Matteo Fischetti. {0,1/2}-Chvatal-Gomory cuts. Math. Programming, Volume 74, p221--235, 1996.
 * - Arie M. C. A. Koster, Adrian Zymolka and Manuel Kutschka. \n
 *   Algorithms to separate {0,1/2}-Chvatal-Gomory cuts. 
 *   Algorithms - ESA 2007: 15th Annual European Symposium, Eilat, Israel, October 8-10, 2007, \n
 *   Proceedings. Lecture Notes in Computer Science, Volume 4698, p. 693--704, 2007.
 * - Arie M. C. A. Koster, Adrian Zymolka and Manuel Kutschka. \n
 *   Algorithms to separate {0,1/2}-Chvatal-Gomory cuts (Extended Version). \n
 *   ZIB Report 07-10, Zuse Institute Berlin, 2007. http://www.zib.de/Publications/Reports/ZR-07-10.pdf
 * - Manuel Kutschka. Algorithmen zur Separierung von {0,1/2}-Schnitten. Diplomarbeit. Technische Universitaet Berlin, 2007.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "string.h"
#include "sepa_zerohalf.h"
#include "scip/buffer.h"
#include "scip/cons_linear.h"
#include "scip/lp.h"
#include "scip/scipdefplugins.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"


#define SEPA_NAME              "zerohalf"
#define SEPA_DESC              "{0,1/2}-cuts separator"
#define SEPA_PRIORITY             -6000
#define SEPA_FREQ                    -1 
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP           TRUE
#define SEPA_DELAY                FALSE

#define DEFAULT_MAXROUNDS             5 /**< maximal number of zerohalf separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        10 /**< maximal number of zerohalf separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          50 /**< maximal number of {0,1/2}-cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     500 /**< maximal number of {0,1/2}-cuts separated per separation round in root node */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */
#define DEFAULT_DECOMPOSEPROBLEM  FALSE /**< should problem be decomposed into subproblems (if possible)? */
#define DEFAULT_MAXDEPTH             -1 /**< separating cuts only if depth <= maxdepth (-1: unlimited) */
#define DEFAULT_MINVIOLATION       0.30 /**< minimal violation of a {0,1/2}-cut to be separated */
#define DEFAULT_FORCECUTSTOLP     FALSE /**< should the cuts be forced to enter the LP? (bypassing SCIPefficacy criteria) */
#define DEFAULT_FORCECUTSTOSEPASTORE FALSE /**< should the cuts be forced to enter SCIP's sepastore?
                                            *   (bypassing SCIPefficicacy criteria, if no other cut is found) */ 
#define DEFAULT_MAXCUTS             100 /**< maximal number of {0,1/2}-cuts determined per separation round 
                                         *  (this includes separated but inefficacious cuts) */
#define DEFAULT_MAXCUTSROOT        1000 /**< maximal number of {0,1/2}-cuts determined per separation round
                                         *   in the root node (this includes separated but inefficacious cuts) */
#define DEFAULT_SUBSCIPOBJECTIVE    'v' /**< auxiliary IP objective function type */
#define DEFAULT_RELAXCONTVARS    FALSE /**< should continuous variables be relaxed by adding variable bounds? */
#define DEFAULT_SCALEFRACCOEFFS   TRUE /**< should rows be scaled to make fractional coefficients integer? */
#define DEFAULT_SUBSCIPSETTINGS    "-" /**< optional settings file of the auxiliary IP (-: none) */
#define DEFAULT_SUBSCIPSOLLIMIT     -1 /**< limits/solutions setting of the auxiliary IP */
#define DEFAULT_SUBSCIPUSEALLSOLS TRUE /**< should all (proper) solutions of the auxiliary IP be used to generate
                                        *   cuts instead of using only the best? */
#define DEFAULT_PPDELTA          0.500 /**< value of delta parameter used in preprocessing method 'd' */
#define DEFAULT_SUBSCIPOBJPEN    0.001 /**< penalty factor used with objective function 'p' of auxiliary IP */

#define DEFAULT_PPMETHODS     "CXGXIM" /**< preprocessing methods and ordering */
#define DEFAULT_SEPAMETHODS       "2g" /**< preprocessing methods and ordering */
#define DEFAULT_MAXNCALLS         -1LL /**< maximal number of calls (-1: unlimited) */
#define DEFAULT_IGNOREPREVIOUSZHCUTS FALSE /**< should zerohalf cuts found in previous callbacks ignored? */
#define DEFAULT_ONLYORIGROWS     FALSE /**< should only original LP rows be considered (i.e. ignore previously added LP rows)? */
#define DEFAULT_USEZHCUTPOOL      TRUE /**< should zerohalf cuts be filtered using a cutpool */

#define DEFAULT_MAXTESTDELTA        10 /**< maximal number of different deltas to try for cmir (-1: unlimited, 0: delta=1) */
#define DEFAULT_TRYNEGSCALING     TRUE /**< should negative values also be tested in scaling for cmir? */

/* cut pool management */
#define ORTHOFUNC                   'e' 
#define MINORTHO                    0.5 

/* SCIPcalcMIR parameters */
#define NNONZOFFSET                 500 
#define BOUNDSWITCH                0.50
#define USEVBDS                    TRUE
#define ALLOWLOCAL                 TRUE
#define FIXINTEGRALRHS             TRUE
#define BOUNDSFORTRANS             NULL
#define BOUNDTYPESFORTRANS         NULL
#define MAXWEIGHTRANGE         10000.00
#define MINFRAC                    0.05
#define MAXFRAC                    1.00

/* SCIPcalcRowIntegralScalar parameters */
#define MAXDNOM                    1000
#define MAXSCALE                 1000.0

/* should variable bounds be used for substituting continuous variables */
#define USEVARBOUNDS               TRUE 


/* --------------------------------------------------------------------------------------------------------------------
 * definitions of enums and some related strings
 * -------------------------------------------------------------------------------------------------------------------- */

/** preprocessing methods, usable within the ppmethods parameter */
enum preprocessingmethods 
   {
      MODGAUSSIANELIMINATION =                'G',
      DELETEZEROROWS =                        'Z',
      DELETEZEROCOLS =                        'z',
      DELETECOLSINGLETONS =                   's',
      ADDTRIVIALCUTS =                        'X',
      DELETEIDENTROWS =                       'I',
      MERGEIDENTCOLS =                        'i',
      DELETELARGESLACKROWS =                  'L',
      DELETESMALLFRACSOLCOLS =                'd',
      DELETEROWSWRTMINSLACK =                 'M',
      PPCOLUMNS =                             'C',
      PPROWS =                                'R'
   };
#if 0 /* currently not used */
typedef enum preprocessingmethods PREPROCESSINGMETHODS;  
#endif

/** separation methods, usable within the sepamethods parameter  */
enum sepamethods 
   {
      STOPIFCUTWASFOUND =                     '!',
      SOLVEAUXSCIP =                          's',
      SOLVEAUXSCIPEXACT =                     'S',
      ENUMHEURNMAX1 =                         'e',
      ENUMHEURNMAX2 =                         'E',
      GAUSSHEUR =                             'g',
      MAX2ODDENTRIESPERROW =                  '2'
   };
#if 0 /* currently not used */
typedef enum sepamethods SEPAMETHODS;
#endif

/** statistics: "origin" of separated cut */
enum cutseparatedby 
   {                                                
      AUXIP, DECOMPOSITION, PPZEROONEROW, HEURISTICSENUM, HEURISTICSGAUSS, AUXGRAPH
   };
typedef enum cutseparatedby CUTSEPARATEDBY;
    


/* --------------------------------------------------------------------------------------------------------------------
 * auxiliary (inline) functions
 * -------------------------------------------------------------------------------------------------------------------- */

#define ISEVEN(scip, value) (SCIPisEQ((scip) , SCIPfloor(scip , (value) / 2) , (value) / 2))      /**< is value even? */
#define ISODD(scip, value) (!(ISEVEN((scip), (value))))                                                    /**< is value odd? */
#define XOR(bool1, bool2) ((bool1) ^ (bool2))
#define DIV(value1, powerof2value) (((unsigned int)(value1)) / ((unsigned int)(powerof2value)))     /**< integer division using a power of 2 as divisor */
#define MOD(value1, powerof2value) (((unsigned int)(value1)) % ((unsigned int)(powerof2value)))     /**< remainder of integer division using a power of 2
                                                                                                       as divisor */
#ifndef BMSmoveMemoryArray                                      
/** moves array at source with size num to ptr */
#define BMSmoveMemoryArray(ptr, source, num)                    \
   {                                                            \
      size_t size__ = (num) * sizeof(*(ptr));                   \
      if( size__ > 0 )                                          \
      {                                                         \
         assert((void*)(ptr) != NULL);                          \
         assert((void*)(source) != NULL);                       \
         memmove((void*)(ptr), (void*)(source), size__);        \
      }                                                         \
   }
#endif

#ifdef  ZEROHALF__PRINT_STATISTICS

#define ZEROHALFstatistics(x)                x                    /**< execute if ZEROHALF__PRINT_STATISTICS is defined */
#define ZEROHALFstatisticsMessage            printf("####   ") ; printf                   /**< print statistics message */
#define ZEROHALFcreateNewTimer(timervar)     SCIP_CALL(SCIPcreateClock(scip, &timervar)) /**< create new timer */
#define ZEROHALFcreateTimer(timervar)        SCIP_CALL(SCIPcreateClock(scip, &timervar)) /**< recreate existing timer */
#define ZEROHALFfreeTimer(timervar)          SCIP_CALL(SCIPfreeClock(scip, &timervar))                /**< free timer */
#define ZEROHALFresetTimer(timervar)         SCIP_CALL(SCIPresetClock(scip, timervar))               /**< reset timer */
#define ZEROHALFstartTimer(timervar)         SCIP_CALL(SCIPstartClock(scip, timervar))               /**< start timer */
#define ZEROHALFstopTimer(timervar)          SCIP_CALL(SCIPstopClock(scip, timervar))                 /**< stop timer */
#define ZEROHALFevalTimer(timervar)          (SCIPgetClockTime(scip, timervar))          /**< evaluate timer (get time) */

#else

#if 0 /* currently not used */
#define ZEROHALFstatistics(x)                /**/                                                          /**< nothing */
#define ZEROHALFstatisticsMessage            while( FALSE ) printf                                         /**< nothing */
#define ZEROHALFcreateNewTimer(timervar)     /**/                                                          /**< nothing */
#define ZEROHALFcreateTimer(timervar)        /**/                                                          /**< nothing */
#define ZEROHALFfreeTimer(timervar)          /**/                                                          /**< nothing */
#define ZEROHALFresetTimer(timervar)         /**/                                                          /**< nothing */
#define ZEROHALFstartTimer(timervar)         /**/                                                          /**< nothing */
#define ZEROHALFstopTimer(timervar)          /**/                                                          /**< nothing */
#define ZEROHALFevalTimer(timervar)          (0.0)                                                         /**< nothing */
#endif

#endif


/* --------------------------------------------------------------------------------------------------------------------
 * symbolic constants for statistical analysis
 * -------------------------------------------------------------------------------------------------------------------- */

/* row / columns */                  
#define IRRELEVANT                          -1        /**< row or column is irrelevant                                  */

/* row */
#define ZERO_ROW                          -100        /**< row has no nonzero entries                                   */
#define IDENT_TO_ROW_WITH_SMALLER_SLACK   -101        /**< row is identical to another row but has a larger slack value */
#define SLACK_GREATER_THAN_MAXSLACK       -102        /**< row has a slack value > maxslack                             */
#define DEFINES_VIOLATED_ZEROHALF_CUT     -103        /**< row defines a violated zerohalf cut                          */
#ifdef WITHDECOMPOSE
#define ROW_IN_SUBPROB_WITHOUT_ODD_RHS    -104        /**< row is part of a subproblem without rows with odd rhs        */
#endif
#define NONEXISTENT_ROW                   -105        /**< row does not exist (lhs is -infinity or rhs is infinity)     */
#define NO_RELEVANT_COLUMNS               -106        /**< row does not contain relevant columns                        */
#define SLACK_GREATER_THAN_MSL_MINUS_SODD -107        /**< row has even rhs and the sum of its slack value and the minimum
                                                         slack value of a odd-rhs-row exceeds maxslack                  */
#define LARGE_COL_EXISTS                  -108        /**< row contains a column which rounding penalty exceeds maxslack
                                                         and the sum of this row's slack and the minimum slack of another
                                                         row with the proper columns exceeds maxslack as well           */

/* column */
#define ZERO_COLUMN                       -200        /**< column has no nonzero entries                                */
#define IDENT_TO_ANOTHER_COLUMN           -201        /**< column is identical to another column                        */
#define ZERO_LP_SOL                       -202        /**< column corresponds to a variable whose LP solution is zero   */
#define LP_SOL_EQUALS_EVEN_LB             -203        /**< column corresponds to a variable whose LP solution equals its even lb */ 
#define LP_SOL_EQUALS_ODD_LB              -204        /**< column corresponds to a variable whose LP solution equals its odd lb  */ 
#define LP_SOL_EQUALS_EVEN_UB             -205        /**< column corresponds to a variable whose LP solution equals its even ub */ 
#define LP_SOL_EQUALS_ODD_UB              -206        /**< column corresponds to a variable whose LP solution equals its odd ub  */ 
#define SINGLETON_COLUMN                  -207        /**< column has only one nonzero entry                            */
#define CONTINUOUS_VARIABLE               -208        /**< column corresponds to a non-integer variable                 */
#define SMALL_FRACSOL_HEUR                -209        /**< column has been omitted (see preprocessColumnsWithSmallFracsol)       */
#define ALL_MATRIX_ROWS_DELETED           -210        /**< all rows (of the current subproblem) have been deleted       */
#ifdef WITHDECOMPOSE
#define COLUMN_IN_SUBPROB_WITHOUT_ODD_RHS -211        /**< column is part of a subproblem without a row with odd rhs value       */
#endif


/* --------------------------------------------------------------------------------------------------------------------
 * bit array data structures and functions
 * -------------------------------------------------------------------------------------------------------------------- */


#define BITARRAYBASETYPE                     unsigned int          /**< base type used for the bitarray data structures */
#define BITARRAYBITMASKTYPE                  BITARRAYBASETYPE

/** size of BITARRAYBASETYPE */
static const unsigned int Zerohalf_bitarraybasetypesize = sizeof(BITARRAYBASETYPE);

/** number of bits per BITARRAYBASETYPE */
static const unsigned int Zerohalf_bitarraybasetypesize_nbits = sizeof(BITARRAYBASETYPE) << 3;


#define BITARRAY                             BITARRAYBASETYPE*

/** get the bit mask where the pos-th bit is set */
#define BITMASK(pos)                         ((unsigned int)(1 << (pos)))

/** set the pos-th bit of var */
#define BITSET(var, pos)                     (var) |= BITMASK(pos)

/** is the pos-th bit of var set? */
#define BITISSET(var, pos)                   (var & BITMASK(pos))

/** set the pos-th bit of bitarray barray */
#define BITARRAYBITSET(barray, pos)          BITSET(barray[DIV((pos),Zerohalf_bitarraybasetypesize_nbits)], \
      MOD(pos,Zerohalf_bitarraybasetypesize_nbits))

/** is the pos-th bit of bitarray barray set? */
#define BITARRAYBITISSET(barray, pos)        BITISSET(barray[DIV(pos,Zerohalf_bitarraybasetypesize_nbits)], \
      MOD(pos,Zerohalf_bitarraybasetypesize_nbits))

/** clear bitarray */
#define BITARRAYCLEAR(barray, barraysize)    BMSclearMemoryArray(barray,barraysize)

/** calculates the number of array elements (w.r.t. the bitarray base type) required to create the bitarray */
#define GETREQUIREDBITARRAYSIZE(nvalstostore)                           \
   ((((unsigned int)(nvalstostore)) % (Zerohalf_bitarraybasetypesize_nbits) == 0) \
      ? (((unsigned int)(nvalstostore)) / (Zerohalf_bitarraybasetypesize_nbits)) \
      : ((((unsigned int)(nvalstostore)) / (Zerohalf_bitarraybasetypesize_nbits)) + 1))

/** get the corresponding array element of a bitarray position */
#define GETBITARRAYINDEX(pos)                DIV((pos),Zerohalf_bitarraybasetypesize_nbits)

/** get the bitmask to mask all bits except the pos-th bit of an array element */
#define GETBITARRAYMASK(pos)                 BITMASK(MOD((pos),Zerohalf_bitarraybasetypesize_nbits))

/** apply operation op for all array elements of bitarray barray1 and barray2 */
#define BITARRAYSFOREACH(barray1, barray2, size, op)                    \
   {                                                                    \
      int idx__;                                                        \
      for( idx__ = 0 ; idx__ < (size) ; ++idx__)                        \
      {                                                                 \
         barray2[idx__] op barray1[idx__];                              \
      }                                                                 \
   }

/** barray2 = barray1 XOR barray2 */
#define BITARRAYSXOR(barray1, barray2, size) BITARRAYSFOREACH(barray1,barray2,size,^=)

/** are barray1 and barray2 equal? */
#define BITARRAYSAREEQUAL(barray1, barray2, size)                       \
   (memcmp((void*)(barray1), (void*)(barray2), (size_t)((size) * (Zerohalf_bitarraybasetypesize))) == 0)

#if 0 /* currently not used */
/** clear the pos-th bit of var */
#define BITCLEAR(var, pos)                   (var) &= ~BITMASK(pos)

/** flip the pos-th bit of var */
#define BITFLIP(var, pos)                    (var) ^= BITMASK(pos)

/** clear the pos-th bit of bitarray barray */
#define BITARRAYBITCLEAR(barray, pos)        BITCLEAR(barray[DIV((pos),Zerohalf_bitarraybasetypesize_nbits)], \
      MOD(pos,Zerohalf_bitarraybasetypesize_nbits))

/** flip the pos-th bit of bitarray barray */
#define BITARRAYBITFLIP(barray, pos)         BITFLIP(barray[DIV((pos),Zerohalf_bitarraybasetypesize_nbits)], \
      MOD(pos,Zerohalf_bitarraybasetypesize_nbits))

/** barray2 = barray1 AND barray2 */
#define BITARRAYSAND(barray1, barray2, size) BITARRAYSFOREACH(barray1,barray2,size,&=)

/** barray2 = barray1 OR barray2 */
#define BITARRAYSOR(barray1, barray2, size)  BITARRAYSFOREACH(barray1,barray2,size,|=)

/** barray2 = NOT barray1 */
#define BITARRAYSNOT(barray1, barray2, size) BITARRAYSFOREACH(barray1,barray2,size,= ~)
#endif


/* --------------------------------------------------------------------------------------------------------------------
 * data structures
 * -------------------------------------------------------------------------------------------------------------------- */


/** parameters */
struct SCIP_SepaData
{
 
   int                   maxrounds;          /**< maximal number of {0,1/2} separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of {0,1/2} separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of {0,1/2} cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of {0,1/2} cuts separated per separation round in root node */
   int                   maxdepth;           /**< separating cuts only if depth <= maxdepth (-1: unlimited) */
   SCIP_Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
   SCIP_Bool             decomposeproblem;   /**< should problem be decomposed into subproblems (if possible)? */  
  
   SCIP_Real             minviolation;       /**< minimal violation of a {0,1/2}-cut to be separated */
   SCIP_Bool             forcecutstolp;      /**< should the cuts be forced to enter the LP? */
   SCIP_Bool             forcecutstosepastore;  /**< should the cuts be forced to enter SCIP's sepastore? */
   int                   maxcuts;            /**< maximal number of {0,1/2}-cuts determined per separation round 
                                                (this includes separated but inefficacious cuts) */
   int                   maxcutsroot;        /**< maximal number of {0,1/2}-cuts determined per separation round
                                                in the root node (this includes separated but inefficacious cuts) */  

   char*                 ppmethods;          /**< preprocessing methods */
   char*                 sepamethods;        /**< separation methods */
   int                   nppmethods;         /**< length of ppmethods string */
   int                   nsepamethods;       /**< length of sepamethods string */
  
   int                   subscipsollimit;    /**< value of auxiliary IP / subscip  "limits/sol" */
   SCIP_Bool             subscipuseallsols;  /**< should all known feasible solution of the auxiliary IP be considered? */
   SCIP_Bool             relaxcontvars;      /**< should continuous vars be relaxed by adding varbounds? */
   SCIP_Bool             scalefraccoeffs;    /**< should fractional coeffs be scaled to become integer? */
   char*                 subscipsettings;    /**< optional settings file for the auxiliary IP / subscip */
   char                  subscipobjective;   /**< type of objective function of the auxiliary IP / subscip */
   SCIP_Real             ppdelta;            /**< value of delta parameter used in preprocessing method 'd' */
   SCIP_Real             subscipobjpen;      /**< penalty factor used with objective function 'p' of auxiliary IP */
   SCIP_Longint          maxncalls;          /**< maximal number of callbacks */
   SCIP_Bool             ignoreprevzhcuts;   /**< should zerohalf cuts found within previous callbacks considered as well? */
   SCIP_Bool             onlyorigrows;       /**< should only original LP rows be considered (i.e. ignore previously added LP rows)? */
   SCIP_Bool             usezhcutpool;       /**< should zerohalf cuts be filtered using a cutpool? */
  
   SCIP_Real             maxslack;           /**< initial: 1.0 - 2.0 * minviolation */
   int                   norigrows;          /**< number of original LP rows */
   int*                  origrows;           /**< set of SCIP_ROW->index of all original LP rows */

   int                   maxnnonz;           /**< maximal number of nonzeros allowed in a zerohalf cut */
   int                   maxtestdelta;       /**< maximal number of different deltas to try for cmir (-1: unlimited, 0: delta=1) */
   SCIP_Bool             trynegscaling;      /**< should negative values also be tested in scaling for cmir? */

   /* statistics */
   int                   totalncutsfound;    /**< total number of separated zerohalf cuts, including inefficious ones */
   int                   totalnsepacuts;     /**< total number of separated zerohalf cuts */
   SCIP_CLOCK**          pptimers;           /**< timers of preprocessing methods */
   SCIP_CLOCK**          sepatimers;         /**< timers of separation algorithms */
   SCIP_CLOCK*           dtimer;             /**< timer of decomposition method */
   int*                  nsepacutsalgo;      /**< number zerohalf cuts separated by a specific separation algorithm,
                                              *   including inefficious cuts */
   int*                  nzerohalfcutsalgo;  /**< number zerohalf cuts separated by a specific separation algorithm */
};


/** sub data of the LP or a sub-LP obtained by problem decomposition */
struct Zerohalf_SubLPData
{
   int*                  rrows;              /**< relevant rows (indices of elements in rows array of ZEROHALF_LPDATA) */
   int                   nrrows;             /**< number of relevant rows */
   SCIP_Real*            rrowsrhs;           /**< rhs value of relevant rows; could also be lhs value of the orig row */
   SCIP_Real*            rrowsslack;         /**< slack value of relevant rows */
   int*                  rcols;              /**< relevant columns (indices of elements in cols array of ZEROHALF_LPDATA) */
   int                   nrcols;             /**< number of relevant columns */
   SCIP_Real*            rcolslbslack;       /**< slack value of lower bound constraint: x*_j - lb_j */
   SCIP_Real*            rcolsubslack;       /**< slack value of upper bound constraint: ub_j - x*_j */
};
typedef struct Zerohalf_SubLPData ZEROHALF_SUBLPDATA;


/** LP data */
struct Zerohalf_LPData
{
   SCIP_VAR**            vars;               /**< LP variables */
   SCIP_ROW**            rows;               /**< LP rows */
   SCIP_COL**            cols;               /**< LP columns */
   int                   nvars;              /**< number of LP variables */
   int                   nrows;              /**< number of LP rows */
   int                   ncols;              /**< number of LP columns */
   int                   nvarbounds;         /**< number of variable bounds (-x_j <= -lb_j, x_j <= ub_j) */
  
   ZEROHALF_SUBLPDATA**  subproblems;        /**< decomposed subproblems (subset of the variables, rows and columns above) */
   int                   nsubproblems;       /**< number of subproblems */

   SCIP_Real*            intscalarsleftrow;  /**< array of scalars that would make left half-rows (-a^Tx <= -lhs) rows integral (0.0 if scalar has not been calculated) */
   SCIP_Real*            intscalarsrightrow; /**< array of scalars that would make right half-rows (a^Tx <= rhs) rows integral (0.0 if scalar has not been calculated) */

   /* row related index sets */
   int*                  subproblemsindexofrow;        /**< is rows index relevant? value <0: not relevant,
                                                          value >=0: index of subproblem containing the row */
   int*                  rrowsindexofleftrow;          /**< maps rows index of lhs <= a^Tx <= rhs to rrows index of -a^Tx <= -lhs */
   int*                  rrowsindexofrightrow;         /**< maps rows index of lhs <= a^Tx <= rhs to rrows index of  a^Tx <=  rhs */

   /* col related index sets */
   int*                  subproblemsindexofcol; /**< is cols index relevant? value <0: not relevant
                                                 *   value >=0: index of subproblem containing the column */
   int*                  rcolsindexofcol;    /**< maps cols index to rcols index */

   int*                  bestlbidxofcol;     /**< maps cols index of a continuous variable to the index of its
                                              *   best lower bound (-2: undetermined, -1: lb, >=0: index of vlb)*/
   int*                  bestubidxofcol;     /**< maps cols index of a continuous variable to the index of its
                                              *   best upper bound (-2: undetermined, -1: ub, >=0: index of vub)*/

   /* statistics */
   int                   ndelvarbounds;      /**< number of deleted variable bounds by basic preprocessing */

};
typedef struct Zerohalf_LPData ZEROHALF_LPDATA;


/** data structure to store data of the auxiliary IP:
 *  (1) minimize violation:
 *      min   z :=           s^T v + x^T y
 *      s.t.       (b (mod 2))^T v         - 2q     = 1   "odd rhs"
 *                 (A (mod 2))^T v -     y -    -2r = 0   "column sum"
 *
 *                 v \\in {0,1}^m
 *                 y \\in {0,1}^n
 *                 r \\in Z^n_+ 
 *                 q \\in Z_+  
 *  (2) minimize (weighted) number of aggregated rows
 *      min   z :=           w^T v
 *      s.t.                 s^T v + x^T y          < 1   "feasibility"
 *                 (b (mod 2))^T v         - 2q     = 1   "odd rhs"
 *                 (A (mod 2))^T v -     y -    -2r = 0   "column sum"
 *
 *                 v \\in {0,1}^m
 *                 y \\in {0,1}^n
 *                 r \\in Z^n_+ 
 *                 q \\in Z_+
 *
 *      with w \\in \\Rset^m
 */
struct Zerohalf_AuxIPData
{
   SCIP*                 subscip;            /**< pointer to (sub)SCIP data structure containing the auxiliary IP */
   int                   m;                  /**< number of rows */
   int                   n;                  /**< number of cols */
  
   SCIP_VAR**            v;                  /**< decision variable: 1 iff row is selected for generating a violated zerohalf cut */
   SCIP_VAR**            y;                  /**< auxiliary variable used for calculating the rounding down penalties in the objective */
   SCIP_VAR**            r;                  /**< auxiliary variable used for modelling mod 2 calculus */
   SCIP_VAR*             q;                  /**< auxiliary variable used for modelling mod 2 calculus */

   SCIP_CONS*            feasipcons;         /**< feasibility constraint */
   SCIP_CONS*            oddrhscons;         /**< odd rhs constraint */
   SCIP_CONS**           columnsumcons;      /**< column sum constraints */

   SCIP_Real             timelimit;          /**< value of "limits/time" of subscip */
   SCIP_Real             memorylimit;        /**< value of "limits/memory" of subscip */
   SCIP_Real             objectivelimit;     /**< value of objective limit of subscip */
   int                   nodelimit;          /**< value of "limits/nodes" of subscip */
};
typedef struct Zerohalf_AuxIPData ZEROHALF_AUXIPDATA;


/** data structure to store data (mod 2) densely*/
struct Zerohalf_Mod2Data
{
   ZEROHALF_SUBLPDATA*   relatedsubproblem;  /**< pointer to corresponding subproblem data structure */
  
   BITARRAY*             rows;               /**< dense mod 2 rows */
   BITARRAY*             rowaggregations;    /**< Zerohalf_Mod2Data->rows index set, storing the actual row aggregations */
   SCIP_Bool*            rhs;                /**< TRUE iff corresponding relatedsubproblem->rrowsrhs is odd */

   SCIP_Real*            slacks;             /**< slack value (equal to corresponding relatedsubproblem->rrowsslack value) */
   SCIP_Real*            fracsol;            /**< LP solution value of variable of SCIP_COL* */
  
   int                   nrows;              /**< number of Zerohalf_Mod2Data->rows */
   int                   rowsbitarraysize;   /**< size (w.r.t. bitarray base type) of Zerohalf_Mod2Data->rows array */
   int                   rowaggregationsbitarraysize;  /** size (w.r.t. bitarray base type) of Zerohalf_Mod2Data->rowaggregations array */
   int                   nvarbounds;         /**< number of variable bounds that are part of Zerohalf_Mod2Data->rows */

   int*                  rowsind;            /**< index set of subset of Zerohalf_Mod2Data->rows */
   int*                  colsind;            /**< index set of subset of relatedsubproblem->rcols */
  
   int                   nrowsind;           /**< number of rowsind elements */
   int                   ncolsind;           /**< number of colsind elements */

   /* statistics */
   int*                  rowstatistics;      /**< stores if and why a row was removed in preprocessing */
   int*                  colstatistics;      /**< stores if and why a column was removed in preprocessing */
};
typedef struct Zerohalf_Mod2Data ZEROHALF_MOD2DATA;


/** data structure to store a violated zerohalf cut and related data */
struct Zerohalf_CutData
{
   ZEROHALF_SUBLPDATA*   relatedsubproblem;  /**< pointer to corresponding subproblem data structure */
   ZEROHALF_MOD2DATA*    relatedmod2data;    /**< pointer to corresponding mod 2 data structure */

   SCIP_ROW*             cut;                /**< pointer to SCIP_ROW storing the cut */

   SCIP_Bool             success;            /**< was SCIPcalcMIR successful? */
   SCIP_Bool             isfeasviolated;     /**< is zerohalf cut violated w.r.t. the feasibility tolerance? */
   SCIP_Bool             islocal;            /**< is zerohalf cut only locally valid? */
  
   SCIP_Real             activity;           /**< activity of the zerohalf cut */
   SCIP_Real             rhs;                /**< rhs of the zerohalf cut */
   SCIP_Real             norm;               /**< norm of the nonzero elements of the zerohalf cut */
   SCIP_Real             efficacy;           /**< efficacy of the zerohalf cut */
   SCIP_Real             violation;          /**< violation of the zerohalf cut */
   int                   nnonz;              /**< number of nonzero coefficients of the zerohalf cut */
  
   /* statistics */  
   int                   nrowsincut;         /**< number of LP rows combined into the zerohalf cut */
   int                   nrrowsincut;        /**< number of preprocessed/aggregated LP rows combined into the zerohalf cut */
   CUTSEPARATEDBY        separatedby;        /**< flag to store the method that has separated the zerohalf cut */
   int                   addedtolp;          /**< was the cut added to SCIP? */

};
typedef struct Zerohalf_CutData ZEROHALF_CUTDATA;
  


/** auxiliary graph node data structure */
struct Zerohalf_AuxGraph_Node;
typedef struct Zerohalf_AuxGraph_Node ZEROHALF_AUXGRAPH_NODE;

struct Zerohalf_AuxGraph_Node
{  
   ZEROHALF_AUXGRAPH_NODE** neighbors;       /**< node adjacency list */
   SCIP_Real*            edgeweights;        /**< weights of outgoing edges */
   int*                  relatedrows;        /**< label mapping outgoing edges to mod 2 rows */
   int                   nneighbors;         /**< number of adjacent nodes */
   SCIP_Real             distance;           /**< actual distant from start node (used by Dijkstra)*/
   ZEROHALF_AUXGRAPH_NODE* previous;         /**< previous node in shortest-path-tree (used by Dijkstra) */
};


/** auxiliary graph data structure */
struct Zerohalf_AuxGraph
{
   ZEROHALF_AUXGRAPH_NODE** nodes;           /**< list of all original nodes */
   ZEROHALF_AUXGRAPH_NODE** nodecopies;      /**< list of all copies of original nodes */
   int                   nnodes;             /**< number of original nodes (equals number of copies) */
};
typedef struct Zerohalf_AuxGraph ZEROHALF_AUXGRAPH;


/* --------------------------------------------------------------------------------------------------------------------
 * local methods: create / free data structure
 * -------------------------------------------------------------------------------------------------------------------- */


/** creates and initializes sub LP data structures */
static
SCIP_RETCODE ZerohalfSubLPDataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_SUBLPDATA**  subproblem          /**< pointer to store pointer to created data structure */
   )
{
   assert(scip != NULL);
   assert(subproblem != NULL);

   SCIP_CALL(SCIPallocMemory(scip, subproblem));
   (*subproblem)->rrows = NULL;
   (*subproblem)->rrowsrhs = NULL;
   (*subproblem)->rrowsslack = NULL;
   (*subproblem)->nrrows = 0;
  
   (*subproblem)->rcols = NULL;
   (*subproblem)->rcolslbslack = NULL;
   (*subproblem)->rcolsubslack = NULL;
   (*subproblem)->nrcols = 0;
  
   return SCIP_OKAY;
}


/** frees sub LP data structures */
static
void ZerohalfSubLPDataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_SUBLPDATA**  subproblem          /**< pointer to pointer of data structure */
   )
{
   assert(scip != NULL);
   assert(subproblem != NULL);
   assert(*subproblem != NULL);

   if( (*subproblem)->rrows != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*subproblem)->rrows));  
   }
   if( (*subproblem)->rrowsrhs != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*subproblem)->rrowsrhs));  
   }
   if( (*subproblem)->rrowsslack != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*subproblem)->rrowsslack));  
   }
   if( (*subproblem)->rcols != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*subproblem)->rcols));  
   }
   if( (*subproblem)->rcolslbslack != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*subproblem)->rcolslbslack)); 
   }
   if( (*subproblem)->rcolsubslack != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*subproblem)->rcolsubslack)); 
   }
   SCIPfreeMemory(scip, subproblem);
   (*subproblem) = NULL;
}


/** creates and initializes LP data structures */
static
SCIP_RETCODE ZerohalfLPDataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_LPDATA**     lpdata              /**< pointer to store pointer to created data structure */
   )
{
   assert(scip != NULL);
   assert(lpdata != NULL);
  
   SCIP_CALL(SCIPallocMemory(scip, lpdata));
   (*lpdata)->vars = NULL;
   (*lpdata)->rows = NULL;
   (*lpdata)->cols = NULL;

   (*lpdata)->subproblems = NULL;
   (*lpdata)->nsubproblems = 0;

   (*lpdata)->intscalarsleftrow = NULL;
   (*lpdata)->intscalarsrightrow = NULL;
   
   (*lpdata)->subproblemsindexofrow = NULL;
   (*lpdata)->rrowsindexofleftrow = NULL;
   (*lpdata)->rrowsindexofrightrow = NULL;
  
   (*lpdata)->subproblemsindexofcol = NULL;
   (*lpdata)->rcolsindexofcol = NULL;  

   (*lpdata)->bestlbidxofcol = NULL;
   (*lpdata)->bestubidxofcol = NULL;

   return SCIP_OKAY;
}


/** frees LP data structures */
static
SCIP_RETCODE ZerohalfLPDataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_LPDATA**     lpdata              /**< pointer to pointer of data structure */
   )
{
   int sp;

   assert(scip != NULL);
   assert(lpdata != NULL);
   assert(*lpdata != NULL);

   /* ! Do not free (*lpdata)->vars, (*lpdata)->rows and (*lpdata)->cols ! */
   assert(((*lpdata)->nsubproblems == 0 && (*lpdata)->subproblems == NULL)
      || ((*lpdata)->nsubproblems > 0 && (*lpdata)->subproblems != NULL));

   if( (*lpdata)->subproblems != NULL )
   {
      for( sp = 0 ; sp < (*lpdata)->nsubproblems ; ++sp )
         if( (*lpdata)->subproblems[sp] != NULL )
            ZerohalfSubLPDataFree(scip, &((*lpdata)->subproblems[sp]));       
      SCIPfreeMemoryArray(scip, &((*lpdata)->subproblems));
   }

   if( (*lpdata)->intscalarsleftrow != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*lpdata)->intscalarsleftrow));
   }
   if( (*lpdata)->intscalarsrightrow != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*lpdata)->intscalarsrightrow));
   }
   if( (*lpdata)->subproblemsindexofrow != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*lpdata)->subproblemsindexofrow));
   }   
   if( (*lpdata)->rrowsindexofleftrow != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*lpdata)->rrowsindexofleftrow));  
   }
   if( (*lpdata)->rrowsindexofrightrow != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*lpdata)->rrowsindexofrightrow));  
   }
   if( (*lpdata)->subproblemsindexofcol != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*lpdata)->subproblemsindexofcol));  
   }
   if( (*lpdata)->rcolsindexofcol != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*lpdata)->rcolsindexofcol));  
   }
   if( (*lpdata)->bestlbidxofcol != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*lpdata)->bestlbidxofcol));  
   }
   if( (*lpdata)->bestubidxofcol != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*lpdata)->bestubidxofcol));  
   }
   SCIPfreeMemory(scip, lpdata);
   (*lpdata) = NULL;
  
   return SCIP_OKAY;
}


/** creates and initializes mod 2 data structures */
static
SCIP_RETCODE ZerohalfMod2DataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_MOD2DATA**   mod2data            /**< pointer to store pointer to created data structure */
   )
{
   assert(scip != NULL);
   assert(mod2data != NULL);
  
   SCIP_CALL(SCIPallocMemory(scip, mod2data));
  
   (*mod2data)->relatedsubproblem = NULL;
  
   (*mod2data)->rows = NULL;
   (*mod2data)->rowaggregations = NULL;
   (*mod2data)->rhs = NULL;

   (*mod2data)->slacks = NULL;
   (*mod2data)->fracsol = NULL;

   (*mod2data)->rowstatistics = NULL;
   (*mod2data)->colstatistics = NULL;  

   (*mod2data)->rowsind = NULL;
   (*mod2data)->colsind = NULL;
  
   return SCIP_OKAY;
}


/** frees  data structures */
static
SCIP_RETCODE ZerohalfMod2DataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_MOD2DATA**   mod2data            /**< pointer to pointer of data structure */
   )
{
   int i;

   assert(scip != NULL);
   assert(mod2data != NULL);
   assert(*mod2data != NULL);

   if( (*mod2data)->rows != NULL )
   {
      for( i = 0 ; i < (*mod2data)->nrows ; ++i)
         if( (*mod2data)->rows[i] != NULL )
         {
            SCIPfreeMemoryArray(scip, &((*mod2data)->rows[i]));
         }
      SCIPfreeMemoryArray(scip, &((*mod2data)->rows));
   }

   if( (*mod2data)->rowaggregations != NULL )
   {
      for( i = 0 ; i < (*mod2data)->nrows ; ++i)
         if( (*mod2data)->rowaggregations[i] != NULL )
         {
            SCIPfreeMemoryArray(scip, &((*mod2data)->rowaggregations[i]));
         }
      SCIPfreeMemoryArray(scip, &((*mod2data)->rowaggregations));
   }

   if( (*mod2data)->rhs != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*mod2data)->rhs));
   }
   if( (*mod2data)->slacks != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*mod2data)->slacks));
   }
   if( (*mod2data)->fracsol != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*mod2data)->fracsol));
   }
   if( (*mod2data)->rowstatistics != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*mod2data)->rowstatistics));
   }
   if( (*mod2data)->colstatistics != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*mod2data)->colstatistics));
   }
   if( (*mod2data)->rowsind != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*mod2data)->rowsind));
   }
   if( (*mod2data)->colsind != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*mod2data)->colsind));
   }
   SCIPfreeMemory(scip, mod2data);
   (*mod2data) = NULL;
  
   return SCIP_OKAY;
}


/** creates and initializes auxiliary IP data structures */
static
SCIP_RETCODE ZerohalfAuxIPDataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_AUXIPDATA**  auxipdata           /**< pointer to store pointer to created data structure */
   )
{  
   assert(scip != NULL);
   assert(auxipdata != NULL);

   SCIP_CALL(SCIPallocMemory(scip, auxipdata));

   (*auxipdata)->subscip = NULL;
  
   (*auxipdata)->v = NULL;
   (*auxipdata)->y = NULL;
   (*auxipdata)->r = NULL;
   (*auxipdata)->q = NULL;

   (*auxipdata)->feasipcons = NULL;    
   (*auxipdata)->oddrhscons = NULL;
   (*auxipdata)->columnsumcons = NULL;
  
   return SCIP_OKAY;
} 


/** frees auxiliary IP data structures */
static
SCIP_RETCODE ZerohalfAuxIPDataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_AUXIPDATA**  auxipdata           /**< pointer to pointer of data structure */
   )
{
   assert(scip != NULL);
   assert(auxipdata != NULL);
   assert(*auxipdata != NULL);

   if( (*auxipdata)->subscip != NULL )
   {
      SCIP_CALL( SCIPfree(&((*auxipdata)->subscip)) );
   }

   if( (*auxipdata)->v != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*auxipdata)->v));
   }
   if( (*auxipdata)->y != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*auxipdata)->y));
   }
   if( (*auxipdata)->r != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*auxipdata)->r));
   }
   if( (*auxipdata)->columnsumcons != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*auxipdata)->columnsumcons));
   }
   SCIPfreeMemory(scip, auxipdata);
   (*auxipdata) = NULL;
  
   return SCIP_OKAY;
}


/** creates and initializes cut data structures */
static
SCIP_RETCODE ZerohalfCutDataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_CUTDATA**    cutdata,            /**< pointer to pointer of data structure */
   ZEROHALF_SUBLPDATA*   relatedsubproblem,  /**< pointer to corresponding subproblem */
   ZEROHALF_MOD2DATA*    relatedmod2data,    /**< pointer to corresponding mod 2 data */
   int                   nrrowsincut,        /**< number of preprocessed / mod 2 rows in cut */
   int                   nrowsincut,         /**< number of original / LP rows in cut */
   CUTSEPARATEDBY        separatedby         /**< flag storing the method that separated this cut */
   )
{
   assert(scip != NULL);
   assert(cutdata != NULL);
   assert(nrrowsincut >= 0);
   assert(nrowsincut >= 0);
  
   SCIP_CALL(SCIPallocMemory(scip, cutdata));

   (*cutdata)->relatedsubproblem = relatedsubproblem;
   (*cutdata)->relatedmod2data = relatedmod2data;

   (*cutdata)->cut = NULL;
  
   (*cutdata)->success = FALSE;
   (*cutdata)->isfeasviolated = FALSE;
   (*cutdata)->islocal = TRUE;

   (*cutdata)->activity = 0.0;
   (*cutdata)->rhs = 0.0;
   (*cutdata)->norm = 0.0;
   (*cutdata)->efficacy = 0.0;
   (*cutdata)->violation = 0.0;
   (*cutdata)->nnonz = 0;
     
   (*cutdata)->nrrowsincut = nrrowsincut;
   (*cutdata)->nrowsincut = nrowsincut;
   (*cutdata)->separatedby = separatedby;
   (*cutdata)->addedtolp = FALSE;

   return SCIP_OKAY;
}


/** frees cut data structures */
static
SCIP_RETCODE ZerohalfCutDataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_CUTDATA**    cutdata             /**< pointer to pointer of data structure */
   )
{
   assert(scip != NULL);
   assert(cutdata != NULL);
   assert(*cutdata != NULL);

   if( (*cutdata)->cut != NULL )
   {
      SCIP_CALL(SCIPreleaseRow(scip, &((*cutdata)->cut)));
   }
   SCIPfreeMemory(scip, cutdata);
   (*cutdata) = NULL;
  
   return SCIP_OKAY;
}


/** creates and initializes auxiliary graph node data structures */
static
SCIP_RETCODE ZerohalfAuxGraphNodeCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_AUXGRAPH_NODE** node             /**< pointer to store pointer to created data structure */
   )
{
   assert(scip != NULL);
   assert(node != NULL);

   SCIP_CALL(SCIPallocMemory(scip, node));

   (*node)->neighbors = NULL;
   (*node)->edgeweights = NULL;
   (*node)->relatedrows = NULL;
   (*node)->nneighbors = 0;

   (*node)->distance = -1.0;
   (*node)->previous = NULL;
  
   return SCIP_OKAY;
}


/** frees auxiliary graph node data structures */
static
SCIP_RETCODE ZerohalfAuxGraphNodeFree(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_AUXGRAPH_NODE** node             /**< pointer to pointer of data structure */
   )
{
   assert(scip != NULL);
   assert(node != NULL);

   if( (*node)->neighbors != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*node)->neighbors));
   }
   if( (*node)->edgeweights != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*node)->edgeweights));
   }
   if( (*node)->relatedrows != NULL )
   {
      SCIPfreeMemoryArray(scip, &((*node)->relatedrows));
   }
   SCIPfreeMemory(scip, node);
   (*node) = NULL;
  
   return SCIP_OKAY;
}


/** creates and initializes auxiliary graph data structures */
static
SCIP_RETCODE ZerohalfAuxGraphCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_AUXGRAPH**   auxgraph            /**< pointer to store pointer to created data structure */
   )
{
   assert(scip != NULL);
   assert(auxgraph != NULL);

   SCIP_CALL(SCIPallocMemory(scip, auxgraph));

   (*auxgraph)->nodes = NULL;
   (*auxgraph)->nodecopies = NULL;

   (*auxgraph)->nnodes = 0;
  
   return SCIP_OKAY;
}


/** frees auxiliary graph data structures */
static
SCIP_RETCODE ZerohalfAuxGraphFree(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_AUXGRAPH**   auxgraph            /**< pointer to pointer of data structure */
   )
{
   int                   n;

   assert(scip != NULL);
   assert(auxgraph != NULL);
     
   if( (*auxgraph)->nodes != NULL )
   {
      assert((*auxgraph)->nnodes > 0);
      for( n = 0; n < (*auxgraph)->nnodes ; ++n )
         if( (*auxgraph)->nodes[n] != NULL )
         {
            SCIP_CALL(ZerohalfAuxGraphNodeFree(scip, &((*auxgraph)->nodes[n])));
         }
      SCIPfreeMemoryArray(scip, (&(*auxgraph)->nodes));
   }
    
   if( (*auxgraph)->nodecopies != NULL )
   {
      assert((*auxgraph)->nnodes > 0);
      for( n = 0; n < (*auxgraph)->nnodes ; ++n )
         if( (*auxgraph)->nodecopies[n] != NULL )
         {
            SCIP_CALL(ZerohalfAuxGraphNodeFree(scip, &((*auxgraph)->nodecopies[n])));
         }
      SCIPfreeMemoryArray(scip, (&(*auxgraph)->nodecopies));
   }
    
   SCIPfreeMemory(scip, auxgraph);
   (*auxgraph) = NULL;
  
   return SCIP_OKAY;
}


/* --------------------------------------------------------------------------------------------------------------------
 * local methods: debug
 * -------------------------------------------------------------------------------------------------------------------- */

#ifdef SCIP_DEBUG
/** returns a string containing the name of the symbolic constant (given as int value) */
static
char* getconstantname(
   char*                 buffer,             /**< string containing the name */
   int                   value               /**< symbolic constant given as int value */
   )
{
   switch( value )
   {
   case IRRELEVANT: 
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "IRRELEVANT"); break;
   case ZERO_ROW: 
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "ZEROROW"); break;
   case IDENT_TO_ROW_WITH_SMALLER_SLACK: 
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "IDENTROW"); break;
   case SLACK_GREATER_THAN_MAXSLACK:
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "SLACK>MAXSLACK"); break;
   case DEFINES_VIOLATED_ZEROHALF_CUT:
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "ISZHCUT"); break;
#ifdef WITHDECOMPOSE
   case ROW_IN_SUBPROB_WITHOUT_ODD_RHS: 
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "NOODDRHSROW"); break;
#endif
   case NONEXISTENT_ROW: 
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "DOESNOTEXIST"); break;
   case NO_RELEVANT_COLUMNS:
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "HASNORELEVANTCOLS"); break;
   case SLACK_GREATER_THAN_MSL_MINUS_SODD: 
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "SLACK+SODD>MAXSLACK"); break;
   case LARGE_COL_EXISTS: 
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "SCOL+SODD>MAXSLACK"); break; 
   case ZERO_COLUMN: 
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "ZEROCOLUMN"); break;
   case IDENT_TO_ANOTHER_COLUMN:
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "IDENTCOLUMN"); break;
   case ZERO_LP_SOL: 
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "PRIMSOL=0"); break;
   case LP_SOL_EQUALS_EVEN_LB:
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "PRIMSOL=EVENLB"); break;
   case LP_SOL_EQUALS_ODD_LB:
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "PRIMSOL=ODDLB"); break;
   case LP_SOL_EQUALS_EVEN_UB:
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "PRIMSOL=EVENUB"); break;
   case LP_SOL_EQUALS_ODD_UB:
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "PRIMSOL=ODDUB"); break;
   case SINGLETON_COLUMN: 
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "SINGLETONCOLUMN"); break;
   case CONTINUOUS_VARIABLE: 
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "CONTCOLUMN"); break;
   case SMALL_FRACSOL_HEUR: 
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "SUMFRACSOL<DELTA"); break;
   case ALL_MATRIX_ROWS_DELETED:
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "NOROWSLEFT"); break;
#ifdef WITHDECOMPOSE
   case COLUMN_IN_SUBPROB_WITHOUT_ODD_RHS:
      (void) SCIPsnprintf(buffer, SCIP_MAXSTRLEN, "NOODDRHSCOL"); break;
#endif
   default: 
      SCIPerrorMessage("parameter <%s> unknown\n", value);
      SCIPABORT(); 
   }

   return buffer;
}
#endif


#ifdef ZEROHALF__PRINT_STATISTICS
/** prints the preprocessing statistics of the basic preprocessing applied while
    reading the LP data and building basic data structures */ 
static
SCIP_RETCODE printPreprocessingStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_LPDATA*      lpdata              /**< data of current LP relaxation */
   )
{
   int                   i;
   int                   nrelevantrows;
   int                   nirrelevantrows;
   int                   nerrors;
   int                   nrelevantcols;
   int                   nirrelevantcols;
   int                   nrows;
   int                   ncols;
   int                   nnonexistingrows;

   assert(lpdata != NULL);

   nrelevantrows    = 0;
   nirrelevantrows  = 0;
   nnonexistingrows = 0;
   nrelevantcols    = 0;
   nirrelevantcols  = 0;
   nerrors          = 0;

   for( i = 0 ; i < lpdata->nrows ; ++i)
   {
      if( lpdata->rrowsindexofleftrow[i] >= 0 )
         nrelevantrows++;
      else
         if( lpdata->rrowsindexofleftrow[i] < -199 )
            nerrors++;
         else
            if( lpdata->rrowsindexofleftrow[i] == NONEXISTENT_ROW )
               nnonexistingrows++;
            else
               nirrelevantrows++;
      if( lpdata->rrowsindexofrightrow[i] >= 0 )
         nrelevantrows++;
      else
         if( lpdata->rrowsindexofrightrow[i] < -199 )
            nerrors++;
         else
            if( lpdata->rrowsindexofrightrow[i] == NONEXISTENT_ROW )
               nnonexistingrows++;
            else
               nirrelevantrows++;
   }
  
   for( i = 0 ; i < lpdata->ncols ; ++i)
   {
      if( lpdata->rcolsindexofcol[i] >= 0 )
         nrelevantcols++;
      else
         if( lpdata->rcolsindexofcol[i] == -1 )
            nirrelevantcols++;
         else
            if( lpdata->rcolsindexofcol[i] > -200 || lpdata->rcolsindexofcol[i] < -299 )
               nerrors++;
            else
               nirrelevantcols++;
   }

   nrows = nrelevantrows + nirrelevantrows - nnonexistingrows;
   ncols = nrelevantcols + nirrelevantcols;
  
   ZEROHALFstatisticsMessage("\n");
   ZEROHALFstatisticsMessage("                | ----- lp data ----- | --- (reductions) -- | --- problem data -- | -lpdata- | -(red.)- | -probd.- | --------\n");
   ZEROHALFstatisticsMessage("                |    nrows |    ncols | ndelrows | ndelcols |   nrrows |   nrcols | nvarbnds | ndlvbnds | nvarbnds |         \n");
   ZEROHALFstatisticsMessage("%15s | %8d | %8d | %8d | %8d | %8d | %8d | %8d | %8d | %8d |\n",
      "READING LPDATA", nrows, ncols, nirrelevantrows - nnonexistingrows, nirrelevantcols, nrelevantrows, nrelevantcols,
      lpdata->nvarbounds, lpdata->ndelvarbounds, lpdata->nvarbounds - lpdata->ndelvarbounds); 

   return SCIP_OKAY;
}
#endif


#ifdef SCIP_DEBUG
/** prints the considered subproblem */
static
void debugPrintSubLpData(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   ZEROHALF_SUBLPDATA*   sublpdata           /**< considered subproblem */
   )
{
   int                   i;
   int                   j;

   assert(scip != NULL);
   assert(lpdata != NULL);
   assert(sublpdata != NULL);

  
   SCIPdebugMessage("\n debugPrintSubLpData:\n\n");
  
   SCIPdebugMessage(" rrows:   (nrrows=%d)\n", sublpdata->nrrows);
   for( i = 0 ; i < sublpdata->nrrows ; ++i)
   {
      SCIPdebugMessage(" %6d:  rrows: %6d  rhs: %6g  slack: %6g  name: %s\n",
         i, sublpdata->rrows[i], sublpdata->rrowsrhs[i], sublpdata->rrowsslack[i],
         SCIProwGetName(lpdata->rows[sublpdata->rrows[i]]));
   }
   SCIPdebugMessage("\n rcols:   (nrcols=%d)\n", sublpdata->nrcols);
   for( j = 0 ; j < sublpdata->nrcols ; ++j)
   {
      SCIPdebugMessage(" %6d:  rcols: %6d  lbslack: %6g  ubslack: %6g\n",
         i, sublpdata->rcols[i], sublpdata->rcolslbslack[i], sublpdata->rcolsubslack[i]);        
   }
}
#endif


#ifdef SCIP_DEBUG
/** prints mod 2 data structures */
static
void debugPrintMod2Data(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered mod 2 data structure */
   SCIP_Bool             printaggregations   /**< should row aggregation bitarrays be printed? */ 
   )
{
   int                   i;
   int                   j;
   int                   k;
  
   assert(scip != NULL);
   assert(lpdata != NULL);
   assert(mod2data != NULL);


   SCIPdebugMessage("\n debugPrintMod2Data:\n\n");

   SCIPdebugMessage(" nrows = %d, nvarbounds = %d, nrcols = %d, nrowsind = %d, ncolsind = %d\n",
      mod2data->nrows, mod2data->nvarbounds, mod2data->relatedsubproblem->nrcols,
      mod2data->nrowsind, mod2data->ncolsind);
   SCIPdebugMessage(" rowsbitarraysize = %d, rowaggregationsbitarraysize = %d\n",
      mod2data->rowsbitarraysize, mod2data->rowaggregationsbitarraysize);
 
  
   SCIPdebugMessage("\n fracsol:\n");
   for( j = 0 ; j < mod2data->relatedsubproblem->nrcols ; ++j )
   {
      for( k = 0 ; k < mod2data->ncolsind ; ++k )
         if( mod2data->colsind[k] == j )
            break;
      SCIPdebugMessage(" rcols[%6d]:  fracsol: %6g  colsind: %6d  name: %s\n", j, mod2data->fracsol[j],
         k < mod2data->ncolsind ? k : -1,
         SCIPvarGetName(SCIPcolGetVar(lpdata->cols[mod2data->relatedsubproblem->rcols[j]])));
   }

   SCIPdebugMessage("\n (A mod 2, b mod 2, [#nonz] (slacks), R, name(rrows), left(-)/right(+):\n");
   if( mod2data->nrowsind == 0 )
   {
      SCIPdebugMessage(" empty\n");
   }
   for( i = 0 ; i < mod2data->nrowsind ; ++i )
   {
      int nnonz = 0;
      SCIPdebugMessage(" ");
      for( j = 0 ; j < mod2data->ncolsind; ++j )
         if( BITARRAYBITISSET(mod2data->rows[mod2data->rowsind[i]], mod2data->colsind[j]) ) /*lint !e701*/
         {
            nnonz++;
            SCIPdebugPrintf("1");        
         }
         else
         {
            SCIPdebugPrintf(".");
         }
      if( mod2data->rhs[mod2data->rowsind[i]] )
      {
         SCIPdebugPrintf("  1");
      }
      else
      {
         SCIPdebugPrintf("  0");
      }
      SCIPdebugPrintf("  [%4d] ", nnonz);    
      SCIPdebugPrintf("(%6g)  ", mod2data->slacks[mod2data->rowsind[i]]);

      if( printaggregations )
      {
         for( j = 0 ; j < mod2data->relatedsubproblem->nrrows; ++j )
            if( BITARRAYBITISSET(mod2data->rowaggregations[mod2data->rowsind[i]], j) ) /*lint !e701*/
            {
               SCIPdebugPrintf("1");
            }
            else
            {
               SCIPdebugPrintf(".");
            }
      }
    
      if( mod2data->rowsind[i] < mod2data->nrows - mod2data->nvarbounds )
      {      
         SCIPdebugPrintf("  %s ", SCIProwGetName(lpdata->rows[mod2data->relatedsubproblem->rrows[mod2data->rowsind[i]]]));
         if( lpdata->rrowsindexofleftrow[mod2data->relatedsubproblem->rrows[mod2data->rowsind[i]]] == mod2data->rowsind[i] )
         {
            SCIPdebugPrintf(" -\n");
         }
         else
         {
            SCIPdebugPrintf(" +\n");
         }
      }
      else
      {
         SCIPdebugPrintf("    varbound(rows[%d])\n", mod2data->rowsind[i]);
      }
   }
   SCIPdebugPrintf("\n");
}
#endif


#ifdef SCIP_DEBUG
/** prints LP data */
static
SCIP_RETCODE debugPrintLPRowsAndCols(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_LPDATA*      lpdata              /**< data of current LP relaxation */
   )
{
   int                   i;
   int                   j;
   char                  temp[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(lpdata != NULL);

   SCIPdebugMessage("\n\nLP rows:\n");
   for( i = 0 ; i < lpdata->nrows ; ++i)
   {
      SCIPdebugMessage("\nrow %d (left): %s[%d,%d] %s:\n", i,
         (lpdata->subproblemsindexofrow[i] == IRRELEVANT)
         || (lpdata->rrowsindexofleftrow[i] < 0) ? "IRRELEVANT" : "RELEVANT",
         lpdata->subproblemsindexofrow[i], lpdata->rrowsindexofleftrow[i],
         lpdata->rrowsindexofleftrow[i] < 0 ? getconstantname(temp, lpdata->rrowsindexofleftrow[i]) : "");
      SCIPdebugMessage("row %d (right): %s[%d,%d] %s:\n", i,
         (lpdata->subproblemsindexofrow[i] == IRRELEVANT)
         || (lpdata->rrowsindexofrightrow[i] < 0) ? "IRRELEVANT" : "RELEVANT",
         lpdata->subproblemsindexofrow[i], lpdata->rrowsindexofrightrow[i],
         lpdata->rrowsindexofrightrow[i] < 0 ? getconstantname(temp, lpdata->rrowsindexofrightrow[i]) : "");
      SCIP_CALL( SCIPprintRow(scip, lpdata->rows[i], NULL) );
   }
  
   SCIPdebugMessage("\n\nLP cols:\n");
   for( j = 0 ; j < lpdata->ncols ; ++j)
   {
      SCIPdebugMessage("\ncol %d: %s[%d,%d] %s:\n", j,
         (lpdata->subproblemsindexofcol[j] == IRRELEVANT)
         || (lpdata->rcolsindexofcol[j] < 0) ? "IRRELEVANT" : "RELEVANT",
         lpdata->subproblemsindexofcol[j], lpdata->rcolsindexofcol[j],
         lpdata->rcolsindexofcol[j] < 0 ? getconstantname(temp, lpdata->rcolsindexofcol[j]) : ""); 
      SCIPdebugMessage("%s = %f\n", 
         SCIPvarGetName(SCIPcolGetVar(lpdata->cols[j])),
         SCIPcolGetPrimsol(lpdata->cols[j]));    
   }

   return SCIP_OKAY;
}
#endif


/* --------------------------------------------------------------------------------------------------------------------
 * local methods: SCIPsortInd comparators
 * -------------------------------------------------------------------------------------------------------------------- */


/** comparator function for sorting an index array non-decreasingly according to a real array */
static
SCIP_DECL_SORTINDCOMP(compRealNonDecreasing)
{
   SCIP_Real* scores;

   scores = (SCIP_Real*) dataptr;

   if( scores[ind1] < scores[ind2] )
      return -1;
   else if( scores[ind1] > scores[ind2] )
      return +1;
   else
      return 0;
}


/** comparator function for sorting an index array non-increasingly according to a real array */
static
SCIP_DECL_SORTINDCOMP(compRealNonIncreasing)
{
   SCIP_Real* scores;
 
   scores = (SCIP_Real*) dataptr;

   if( scores[ind1] < scores[ind2] )
      return +1;
   else if( scores[ind1] > scores[ind2] )
      return -1;
   else
      return 0;
}


/* --------------------------------------------------------------------------------------------------------------------
 * local methods: LP data
 * -------------------------------------------------------------------------------------------------------------------- */


/** searches for relevant columns, i.e., columns that cannot be deleted because of basic preprocessing methods */
static
SCIP_RETCODE getRelevantColumns(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_LPDATA*      lpdata              /**< data of current LP relaxation */
   )
{
   SCIP_VAR*             var;
   SCIP_COL*             col;
   SCIP_Real             primsol;
   SCIP_Real             lb;
   SCIP_Real             ub;
   SCIP_Real             lbslack;
   SCIP_Real             ubslack;
   ZEROHALF_SUBLPDATA*   problem;
   int                   j;
#ifdef ZEROHALF__PRINT_STATISTICS       
   int                   tempnvarbnds;
#endif
   int                   nsubproblems;
  
   assert(scip != NULL);
   assert(lpdata != NULL);
   assert(lpdata->cols != NULL);
   assert(lpdata->ncols > 0);
   assert(lpdata->nrows > 0);
   assert(lpdata->subproblems == NULL);
   assert(lpdata->nsubproblems == 0);
   assert(lpdata->subproblemsindexofrow == NULL);
   assert(lpdata->rrowsindexofleftrow == NULL);
   assert(lpdata->rrowsindexofrightrow == NULL);
   assert(lpdata->subproblemsindexofcol == NULL);
   assert(lpdata->rcolsindexofcol == NULL);
   assert(lpdata->bestlbidxofcol == NULL);
   assert(lpdata->bestubidxofcol == NULL);
  
   nsubproblems = 1;

   /* allocate temporary memory for column data structures */
   SCIP_CALL(SCIPallocMemoryArray(scip, &(lpdata->subproblems), nsubproblems)); /* create one "sub"problem */
   SCIP_CALL(SCIPallocMemoryArray(scip, &(lpdata->subproblemsindexofcol), lpdata->ncols));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(lpdata->rcolsindexofcol), lpdata->ncols));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(lpdata->bestlbidxofcol), lpdata->ncols));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(lpdata->bestubidxofcol), lpdata->ncols));

   SCIP_CALL(ZerohalfSubLPDataCreate(scip, &problem));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(problem->rcols), lpdata->ncols));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(problem->rcolslbslack), lpdata->ncols));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(problem->rcolsubslack), lpdata->ncols));
  
   /* initialize data */
   BMSclearMemoryArray(lpdata->subproblemsindexofcol, lpdata->ncols);
   BMSclearMemoryArray(lpdata->rcolsindexofcol, lpdata->ncols);

   lpdata->nsubproblems = 1;
   lpdata->subproblems[0] = problem;
   lpdata->nvarbounds = 0;
   lpdata->ndelvarbounds = 0;
  
   /* check all cols */
   for( j = 0 ; j < lpdata->ncols ; ++j)
   {  
      /* initialize best lb and best ub (-2: undetermined)*/
      lpdata->bestlbidxofcol[j] = -2;
      lpdata->bestubidxofcol[j] = -2;
 
      col = lpdata->cols[j];    
      var = SCIPcolGetVar(col);

      /* check if vartype is continuous */
      if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
      {                   
         primsol = SCIPcolGetPrimsol(col);
         if( SCIPisFeasZero(scip, primsol) )
            primsol = 0.0;

         assert(SCIPisFeasEQ(scip, SCIPgetVarSol(scip, var), primsol));      
      
         lb = SCIPcolGetLb(col);
         ub = SCIPcolGetUb(col);
         lbslack = primsol - lb;
         ubslack = ub - primsol;
         if( SCIPisFeasZero(scip, lbslack) ) 
            lbslack = 0.0;
         if( SCIPisFeasZero(scip, ubslack) ) 
            ubslack = 0.0;
         assert(SCIPisLE(scip, lb, ub));
         assert(!SCIPisNegative(scip, lbslack));
         assert(!SCIPisNegative(scip, ubslack));
      
#ifdef ZEROHALF__PRINT_STATISTICS       
         tempnvarbnds = 0;
         if( !SCIPisInfinity(scip, (-1) * lb) )
            tempnvarbnds++;
         if( !SCIPisInfinity(scip, ub) )
            tempnvarbnds++;
         lpdata->nvarbounds += tempnvarbnds;
#endif      

         if( SCIPisNegative(scip, lb) )
         {
            /* column is declared to be irrelevant because its lower bound is negative, 
             * variable would have to be shifted, complemented or decomposed. 
             */
            /**@todo consider general integers with negative lower bounds and transform to positive representation
             * and propagate through corresponding rows. In the current version, redundant inequalities might be 
             * considered as cut candidates and valid cuts might be missed, but no wrong cuts should be produced 
             * (due to SCIPcalcMIR) this leads to performance deterioration in the (rare) case of general integers 
             * with negative bounds.
             */
            lpdata->subproblemsindexofcol[j] = IRRELEVANT;
            lpdata->rcolsindexofcol[j] = IRRELEVANT;      
         }
         else if( !SCIPisZero(scip, primsol) )
         {
            if( SCIPisZero(scip, lbslack) )
            {
#ifdef ZEROHALF__PRINT_STATISTICS       
               lpdata->ndelvarbounds += tempnvarbnds;
#endif          
               lpdata->subproblemsindexofcol[j] = IRRELEVANT;
               if( ISODD(scip, lb) )
                  lpdata->rcolsindexofcol[j] = LP_SOL_EQUALS_ODD_LB;
               else
                  lpdata->rcolsindexofcol[j] = LP_SOL_EQUALS_EVEN_LB;
            }
            else
            {
               if( SCIPisZero(scip, ubslack) )
               {
#ifdef ZEROHALF__PRINT_STATISTICS       
                  lpdata->ndelvarbounds += tempnvarbnds;
#endif            
                  lpdata->subproblemsindexofcol[j] = IRRELEVANT;
                  if( ISODD(scip, ub) )
                     lpdata->rcolsindexofcol[j] = LP_SOL_EQUALS_ODD_UB;
                  else
                     lpdata->rcolsindexofcol[j] = LP_SOL_EQUALS_EVEN_UB;
               }
               else
               {
                  /* relevant col was found */
                  problem->rcols[problem->nrcols] = j;
                  problem->rcolslbslack[problem->nrcols] = lbslack;
                  problem->rcolsubslack[problem->nrcols] = ubslack;
                  lpdata->subproblemsindexofcol[j] = 0;
                  lpdata->rcolsindexofcol[j] = problem->nrcols;
                  problem->nrcols++;
               }           
            }
         }
         else
         {
#ifdef ZEROHALF__PRINT_STATISTICS       
            lpdata->ndelvarbounds += tempnvarbnds;
#endif        
            lpdata->subproblemsindexofcol[j] = IRRELEVANT;
            lpdata->rcolsindexofcol[j] = ZERO_LP_SOL;
         }
      }
      else
      {
         /* column is irrelevant because vartype is continuous (is handled in getRelevantRows())*/
         lpdata->subproblemsindexofcol[j] = IRRELEVANT;
         lpdata->rcolsindexofcol[j] = CONTINUOUS_VARIABLE;
      }
   }
  
   return SCIP_OKAY;
}


/** finds closest lower bound of col and stores it within lpdata;
 *  the bound can be the lower bound or the best variable lower bound with nonnegative column variable 
 */
static
void findClosestLb(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */   
   SCIP_COL*             col,                /**< column to get closest lower bound */
   SCIP_Real*            bestlbsol,          /**< pointer to store value of closest lower bound */
   int*                  bestlbtype,         /**< pointer to store type of closest lower bound */
   SCIP_VAR**            bestzvlb,           /**< pointer to store variable z in closest variable lower bound b*z + d */
   SCIP_Real*            bestbvlb,           /**< pointer to store coefficient b in closest variable lower bound b*z + d */
   SCIP_Real*            bestdvlb            /**< pointer to store constant d in closest variable lower bound b*z + d */
  
   )
{
   SCIP_VAR* var;
   SCIP_VAR** zvlb;
   SCIP_Real* bvlb;
   SCIP_Real* dvlb;
   int collppos;
   int nvlb;
   int j;

   assert(lpdata != NULL);
   assert(bestlbsol != NULL);
   assert(bestlbtype != NULL);
   assert(bestzvlb != NULL);
   assert(bestbvlb != NULL);
   assert(bestdvlb != NULL);


   collppos = SCIPcolGetLPPos(col);
   var = SCIPcolGetVar(col);
   *bestlbsol = SCIPcolGetLb(col);
   *bestlbtype = lpdata->bestlbidxofcol[collppos];
   *bestzvlb = NULL;
   *bestbvlb = 0.0;
   *bestdvlb = 0.0;

   if( *bestlbtype == -1 )
      return;
   
   if( USEVARBOUNDS ) /*lint !e774 !e506*/
   {
      nvlb = SCIPvarGetNVlbs(var);
      zvlb = SCIPvarGetVlbVars(var);
      bvlb = SCIPvarGetVlbCoefs(var);
      dvlb = SCIPvarGetVlbConstants(var);

      assert(zvlb != NULL || nvlb == 0);
      assert(bvlb != NULL || nvlb == 0);
      assert(dvlb != NULL || nvlb == 0);
   }

   if( *bestlbtype == -2 )
   {
      if( USEVARBOUNDS ) /*lint !e774 !e506*/
      {
         
         /* search for lb or vlb with maximal bound value */
         for( j = 0; j < nvlb; j++ )
         {
            assert(zvlb != NULL);
            assert(bvlb != NULL);
            assert(dvlb != NULL);
            assert(SCIPvarGetType(zvlb[j]) != SCIP_VARTYPE_CONTINUOUS);
         
            /* use only vlb with nonnegative variable z that are column variables and present in the current LP */
            if( SCIPvarGetStatus(zvlb[j]) == SCIP_VARSTATUS_COLUMN && SCIPcolIsInLP(SCIPvarGetCol(zvlb[j])) &&
               !SCIPisNegative(scip, SCIPcolGetLb(SCIPvarGetCol(zvlb[j]))) )
            {
               SCIP_Real vlbsol;
            
               vlbsol = bvlb[j] * SCIPcolGetPrimsol(SCIPvarGetCol(zvlb[j])) + dvlb[j];
               if( vlbsol > *bestlbsol )
               {
                  *bestlbsol = vlbsol;
                  *bestlbtype = j;
               }
            }
         }
      }    

      /* if no better var bound could be found, set type to the fixed bound (-1) */
      if( *bestlbtype == -2 )
         *bestlbtype = -1;      
      
      /* store best bound for substitution */
      lpdata->bestlbidxofcol[collppos] = *bestlbtype;
   }
   assert(lpdata->bestlbidxofcol[collppos] > -2);

   if( *bestlbtype >= 0 )
   {
      assert(USEVARBOUNDS); /*lint !e774 !e506*/
      assert(*bestlbtype < nvlb);
      assert(zvlb != NULL);
      assert(bvlb != NULL);
      assert(dvlb != NULL);
      *bestzvlb = zvlb[*bestlbtype];
      *bestbvlb = bvlb[*bestlbtype];
      *bestdvlb = dvlb[*bestlbtype];      
   }  
}


/** finds closest upper bound of col and stores it within lpdata;
 *  the bound can be the upper bound or the best variable upper bound with nonnegative column variable 
 */
static
void findClosestUb(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */   
   SCIP_COL*             col,                /**< column to get closest upper bound */
   SCIP_Real*            bestubsol,          /**< pointer to store value of closest upper bound */
   int*                  bestubtype,         /**< pointer to store type of closest upper bound */
   SCIP_VAR**            bestzvub,           /**< pointer to store variable z in closest variable upper bound b*z + d */
   SCIP_Real*            bestbvub,           /**< pointer to store coefficient b in closest variable upper bound b*z + d */
   SCIP_Real*            bestdvub            /**< pointer to store constant d in closest variable upper bound b*z + d */

   )
{
   SCIP_VAR* var;
   SCIP_VAR** zvub;
   SCIP_Real* bvub;
   SCIP_Real* dvub;
   int collppos;
   int nvub;
   int j;

   assert(lpdata != NULL);
   assert(bestubsol != NULL);
   assert(bestubtype != NULL);
   assert(bestzvub != NULL);
   assert(bestbvub != NULL);
   assert(bestdvub != NULL);

   collppos = SCIPcolGetLPPos(col);
   var = SCIPcolGetVar(col);
   *bestubsol = SCIPcolGetUb(col);
   *bestubtype = lpdata->bestubidxofcol[collppos];
   *bestzvub = NULL;
   *bestbvub = 0.0;
   *bestdvub = 0.0;

   if( *bestubtype == -1 )
      return;
   
   if( USEVARBOUNDS ) /*lint !e774 !e506*/
   {
      nvub = SCIPvarGetNVubs(var);
      zvub = SCIPvarGetVubVars(var);
      bvub = SCIPvarGetVubCoefs(var);
      dvub = SCIPvarGetVubConstants(var);

      assert(zvub != NULL || nvub == 0);
      assert(bvub != NULL || nvub == 0);
      assert(dvub != NULL || nvub == 0);
   }

   if( *bestubtype == -2 )
   {
      if( USEVARBOUNDS ) /*lint !e774 !e506*/
      {
         /* search for ub or vub with maximal bound value */
         for( j = 0; j < nvub; j++ )
         {
            assert(zvub != NULL);
            assert(bvub != NULL);
            assert(dvub != NULL);
            assert(SCIPvarGetType(zvub[j]) != SCIP_VARTYPE_CONTINUOUS);
         
            /* use only vub with nonnegative variable z that are column variables and present in the current LP */
            if( SCIPvarGetStatus(zvub[j]) == SCIP_VARSTATUS_COLUMN && SCIPcolIsInLP(SCIPvarGetCol(zvub[j])) &&
               !SCIPisNegative(scip, SCIPcolGetUb(SCIPvarGetCol(zvub[j]))) )
            {
               SCIP_Real vubsol;
            
               vubsol = bvub[j] * SCIPcolGetPrimsol(SCIPvarGetCol(zvub[j])) + dvub[j];
               if( vubsol < *bestubsol )
               {
                  *bestubsol = vubsol;
                  *bestubtype = j;
               }
            }
         }
      }
      /* if no better var bound could be found, set type to the fixed bound (-1) */
      if( *bestubtype == -2 )
         *bestubtype = -1;    

      /* store best bound for substitution */
      lpdata->bestubidxofcol[collppos] = *bestubtype;
   }
   assert(lpdata->bestubidxofcol[collppos] > -2);

   if( *bestubtype >= 0 )
   {
      assert(USEVARBOUNDS); /*lint !e774 !e506*/
      assert(*bestubtype < nvub);
      assert(zvub != NULL);
      assert(bvub != NULL);
      assert(dvub != NULL);
      *bestzvub = zvub[*bestubtype];
      *bestbvub = bvub[*bestubtype];
      *bestdvub = dvub[*bestubtype];      
   }  
}


/** searches for relevant rows, i.e., rows containing relevant columns that cannot be deleted because of basic 
 *  preprocessing methods 
 */
static
SCIP_RETCODE getRelevantRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   ZEROHALF_LPDATA*      lpdata              /**< data of current LP relaxation */
   )
{
   SCIP_COL**            colscurrentrow;
   SCIP_Real*            valscurrentrow;
   SCIP_Real*            densecoeffscurrentleftrow;
   SCIP_Real*            densecoeffscurrentrightrow;
   SCIP_ROW*             row;
   SCIP_VAR*             var;
   SCIP_VAR*             bestzvbnd;
   SCIP_Real             bestbndsol;
   SCIP_Real             bestbvbnd;
   SCIP_Real             bestdvbnd;
   SCIP_Real             intscalarleftrow;
   SCIP_Real             intscalarrightrow;
   SCIP_Real             act;
   SCIP_Real             cst;
   SCIP_Real             lhs;
   SCIP_Real             lhsslack;
   SCIP_Real             maxslack;
   SCIP_Real             rhs;
   SCIP_Real             rhsslack;   
   int                   bestbndtype;
   int                   nnonzcurrentrow;
   int                   c;
   int                   r;
   int                   k;
   SCIP_Bool             lhsslackislessequalmaxslack;
   SCIP_Bool             rhsslackislessequalmaxslack;
   SCIP_Bool             lhsisinfinity;
   SCIP_Bool             rhsisinfinity;
   SCIP_Bool             lhsiseven;
   SCIP_Bool             rhsiseven;
   SCIP_Bool             success;
   ZEROHALF_SUBLPDATA*   problem;

   int                   collppos;
   SCIP_Bool             rowisrelevant;
  
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(lpdata != NULL);
   assert(lpdata->rows != NULL);
   assert(lpdata->nrows > 0);
   assert(lpdata->subproblemsindexofcol != NULL);
   assert(lpdata->rcolsindexofcol != NULL);
   assert(lpdata->subproblems != NULL);
   assert(lpdata->nsubproblems > 0);

   assert(lpdata->subproblemsindexofrow == NULL);
   assert(lpdata->rrowsindexofleftrow == NULL);
   assert(lpdata->rrowsindexofrightrow == NULL);

   k = 0;
   problem = lpdata->subproblems[k];

   assert(k >= 0);
   assert(k <= lpdata->nsubproblems);
   assert(problem != NULL);
   assert(problem->rcols != NULL);
   assert(problem->nrcols > 0);
   assert(problem->rcolslbslack != NULL);
   assert(problem->rcolsubslack != NULL);
  
   assert(problem->rrows == NULL);
   assert(problem->rrowsrhs == NULL);
   assert(problem->rrowsslack == NULL);

   /* allocate temporary memory for row data structures */
   SCIP_CALL(SCIPallocMemoryArray(scip, &(lpdata->subproblemsindexofrow), lpdata->nrows));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(lpdata->rrowsindexofleftrow), lpdata->nrows));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(lpdata->rrowsindexofrightrow), lpdata->nrows));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(problem->rrows), 2 * lpdata->nrows));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(problem->rrowsrhs), 2 * lpdata->nrows));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(problem->rrowsslack), 2 * lpdata->nrows));
 
   /* allocate temporary memory */
   SCIP_CALL(SCIPallocBufferArray(scip, &densecoeffscurrentleftrow, lpdata->ncols));
   SCIP_CALL(SCIPallocBufferArray(scip, &densecoeffscurrentrightrow, lpdata->ncols));
   
   /* initialize arrays */
   BMSclearMemoryArray(lpdata->subproblemsindexofrow, lpdata->nrows);
   BMSclearMemoryArray(lpdata->rrowsindexofleftrow, lpdata->nrows);
   BMSclearMemoryArray(lpdata->rrowsindexofrightrow, lpdata->nrows);
   BMSclearMemoryArray(densecoeffscurrentleftrow, lpdata->ncols);
   BMSclearMemoryArray(densecoeffscurrentrightrow, lpdata->ncols);
   
   maxslack = sepadata->maxslack;
   problem->nrrows = 0;
   for( r = 0 ; r < lpdata->nrows ; ++r)
   {     
      row = lpdata->rows[r];

      if( sepadata->ignoreprevzhcuts )
      {
         /* ignore rows whose names start with "zerohalf" */
         const char* rowname = SCIProwGetName(row);
       
         if( strlen(rowname) > 8 )
            if(rowname[0] == 'z'
               && rowname[1] == 'e'
               && rowname[2] == 'r'
               && rowname[3] == 'o'
               && rowname[4] == 'h'
               && rowname[5] == 'a'
               && rowname[6] == 'l'
               && rowname[7] == 'f' )
            {
               lpdata->subproblemsindexofrow[r] = IRRELEVANT;
               lpdata->rrowsindexofleftrow[r] = NONEXISTENT_ROW;
               lpdata->rrowsindexofrightrow[r] = NONEXISTENT_ROW;
               continue;
            }
      }

      /* check if current row is an original LP row (i.e., was not added) if necessary */
      if( sepadata->onlyorigrows )
      {
         int left;
         int center;
         int right;
         int rowindex;
          
         assert(sepadata->origrows != NULL);
         assert(sepadata->norigrows > 0);

         left = 0;
         center = 1;
         right = sepadata->norigrows - 1;
         rowindex = SCIProwGetIndex(row);
         while( left <= right && center > -1 )
         {
            center = left + ((right - left) / 2);
            if( sepadata->origrows[center] == rowindex )
               center = -1;
            if( sepadata->origrows[center] > rowindex )
               right = center - 1;
            else
               left = center + 1;
         }
         if( center > -1 )
         {
            lpdata->subproblemsindexofrow[r] = IRRELEVANT;
            lpdata->rrowsindexofleftrow[r] = NONEXISTENT_ROW;
            lpdata->rrowsindexofrightrow[r] = NONEXISTENT_ROW;
            continue;
         }
      }    
      
      /* get row data */
      colscurrentrow = SCIProwGetCols(row);
      nnonzcurrentrow = SCIProwGetNLPNonz(row);
      valscurrentrow = SCIProwGetVals(row);

      /* clear dense coeffs arrays */
      BMSclearMemoryArray(densecoeffscurrentleftrow, lpdata->ncols);
      BMSclearMemoryArray(densecoeffscurrentrightrow, lpdata->ncols);
      
      /* calculate dense coeffs arrays */
      for( c = 0; c < nnonzcurrentrow; ++c)
      {
         collppos = SCIPcolGetLPPos(colscurrentrow[c]);
	 assert(0 <= collppos && collppos < lpdata->ncols);

         var = SCIPcolGetVar(colscurrentrow[c]);

         /* check if row contains a continuous variable */
         if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         { 
            bestzvbnd = NULL;
            bestbndsol = 0.0;
            bestbvbnd = 0.0;
            bestdvbnd = 0.0;
            bestbndtype = -2;
          
            /* Consider rhs of row and relax continuous variables by substituting for:
             * - a_j > 0: x_j = lb  or  x_j = b*z + d with variable lower bound b*z + d with column var z >= 0
             * - a_j < 0: x_j = ub  or  x_j = b*z + d with variable upper bound b*z + d with column var z >= 0
             * and 
             * consider lhs of row and relax continuous variables by substituting for:
             * - a_j < 0: x_j = lb  or  x_j = b*z + d with variable lower bound b*z + d with column var z >= 0
             * - a_j > 0: x_j = ub  or  x_j = b*z + d with variable upper bound b*z + d with column var z >= 0
             */

            findClosestLb(scip, lpdata, colscurrentrow[c],
               &bestbndsol, &bestbndtype, &bestzvbnd, &bestbvbnd, &bestdvbnd );
            assert( bestbndtype > -2 && lpdata->bestlbidxofcol[collppos] == bestbndtype);
          
            if( bestbndtype > -1 )
            {
               int zlppos;

               zlppos = SCIPcolGetLPPos(SCIPvarGetCol(bestzvbnd));
               assert(0 <= zlppos && zlppos < lpdata->ncols);

               if( valscurrentrow[c] > 0 )
                  densecoeffscurrentrightrow[zlppos] += valscurrentrow[c] * bestbvbnd;
               else
                  densecoeffscurrentleftrow[zlppos] -= valscurrentrow[c] * bestbvbnd;
            }

            findClosestUb(scip, lpdata, colscurrentrow[c],
               &bestbndsol, &bestbndtype, &bestzvbnd, &bestbvbnd, &bestdvbnd );
            assert(bestbndtype > -2 && lpdata->bestubidxofcol[collppos] == bestbndtype);

            if( bestbndtype > -1 )
            {
               int zlppos;

               zlppos = SCIPcolGetLPPos(SCIPvarGetCol(bestzvbnd));
               assert(0 <= zlppos && zlppos < lpdata->ncols);

               if( valscurrentrow[c] > 0 )
                  densecoeffscurrentleftrow[zlppos] -= valscurrentrow[c] * bestbvbnd;
               else
                  densecoeffscurrentrightrow[zlppos] += valscurrentrow[c] * bestbvbnd;
            }          
         }
         else
         {
            densecoeffscurrentleftrow[collppos]  -= valscurrentrow[c];
            densecoeffscurrentrightrow[collppos] += valscurrentrow[c];          
         }
      }
     
      /* calculate scalar that would make (left|right) row coefficients integral; 
       * try to avoid unnecessary or expensive scaling calls 
       */ 
      intscalarleftrow = 1.0;
      intscalarrightrow = 1.0;

      if( sepadata->scalefraccoeffs && (!SCIPisIntegral(scip, SCIPgetRowMinCoef(scip, row)) 
            || !SCIPisIntegral(scip, SCIPgetRowMaxCoef(scip, row)) || !SCIPisIntegral(scip, SCIProwGetSumNorm(row))) )
      {
         SCIP_CALL( SCIPcalcIntegralScalar(densecoeffscurrentleftrow, lpdata->ncols,
               -SCIPepsilon(scip), SCIPepsilon(scip), (SCIP_Longint) MAXDNOM, MAXSCALE, &intscalarleftrow, &success) );
         if( !success )
         {
            lpdata->rrowsindexofleftrow[r] = NONEXISTENT_ROW;
         }
        
         SCIP_CALL( SCIPcalcIntegralScalar(densecoeffscurrentrightrow, lpdata->ncols,
               -SCIPepsilon(scip), SCIPepsilon(scip), (SCIP_Longint) MAXDNOM, MAXSCALE, &intscalarrightrow, &success) );
         if( !success )
         {
            lpdata->rrowsindexofrightrow[r] = NONEXISTENT_ROW;
         }

         if ( lpdata->rrowsindexofleftrow[r] == NONEXISTENT_ROW
            && lpdata->rrowsindexofrightrow[r] == NONEXISTENT_ROW )
         {
            lpdata->subproblemsindexofrow[r] = IRRELEVANT;
            continue;
         }
      }

      lpdata->intscalarsleftrow[r]  = intscalarleftrow;
      lpdata->intscalarsrightrow[r] = intscalarrightrow;

      /* calculate lhs/rhs & slacks */
      act = SCIPgetRowLPActivity(scip, row);
      lhs = SCIProwGetLhs(row);
      rhs = SCIProwGetRhs(row);    
      cst = SCIProwGetConstant(row);
    
      lhsisinfinity = SCIPisInfinity(scip, -lhs);
      rhsisinfinity = SCIPisInfinity(scip, rhs);

      lhsslack = SCIPisFeasZero(scip, act - lhs) ? 0.0 : act - lhs;
      rhsslack = SCIPisFeasZero(scip, rhs - act) ? 0.0 : rhs - act;
    
      lhs = (lhs - cst) * intscalarleftrow;
      rhs = (rhs - cst) * intscalarrightrow;

      lhsisinfinity = lhsisinfinity || SCIPisInfinity(scip, -lhs);
      rhsisinfinity = rhsisinfinity || SCIPisInfinity(scip, rhs);

      lhsslack = lhsslack * intscalarleftrow;
      rhsslack = rhsslack * intscalarrightrow;

      /* check if the slack value of the row is small enough */
      if( (!lhsisinfinity && SCIPisLE(scip, lhsslack, maxslack))
         || (!rhsisinfinity && SCIPisLE(scip, rhsslack, maxslack)) )
      { 
         colscurrentrow = SCIProwGetCols(row);
         nnonzcurrentrow = SCIProwGetNLPNonz(row);
         valscurrentrow = SCIProwGetVals(row);

         lhsiseven = ISEVEN(scip, lhs);
         rhsiseven = ISEVEN(scip, rhs);

         rowisrelevant = FALSE;
         for( c = 0 ; c < nnonzcurrentrow ; ++c )
         {
            collppos = SCIPcolGetLPPos(colscurrentrow[c]);
            var = SCIPcolGetVar(colscurrentrow[c]);

            /* check if row contains a column with primsol = odd bound and update lhs/rhs parity */
            if( lpdata->rcolsindexofcol[collppos] == LP_SOL_EQUALS_ODD_LB
               || lpdata->rcolsindexofcol[collppos] == LP_SOL_EQUALS_ODD_UB )
            {
               assert(SCIPvarGetType(SCIPcolGetVar(colscurrentrow[c])) != SCIP_VARTYPE_CONTINUOUS);
               lhsiseven = !lhsiseven;
               rhsiseven = !rhsiseven;
            }

            /* check if row contains a continuous variable */
            if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
            {   
               if( !rhsisinfinity )
               {
                  if( valscurrentrow[c] * intscalarrightrow > 0 )
                  {
                     findClosestLb(scip, lpdata, colscurrentrow[c], &bestbndsol, &bestbndtype, &bestzvbnd, &bestbvbnd, &bestdvbnd );
                     assert(bestbndtype > -2 && lpdata->bestlbidxofcol[collppos] == bestbndtype);
                  }
                  else
                  {                  
                     assert(valscurrentrow[c] * intscalarrightrow < 0);
                     findClosestUb(scip, lpdata, colscurrentrow[c], &bestbndsol, &bestbndtype, &bestzvbnd, &bestbvbnd, &bestdvbnd );
                     assert(bestbndtype > -2 && lpdata->bestubidxofcol[collppos] == bestbndtype);
                  }
                  assert(bestbndtype == -1 || bestzvbnd != NULL);
                  
                  if( SCIPisInfinity(scip, -bestbndsol) || SCIPisInfinity(scip, bestbndsol) )
                     rhsisinfinity = TRUE;
                  else
                  {
                     /**@todo check whether REALABS is really correct */
                     if ( bestbndtype == -1 )
                        rhs -= intscalarrightrow * REALABS(valscurrentrow[c]) * bestbndsol;
                     else
                        rhs -= intscalarrightrow * REALABS(valscurrentrow[c]) * bestdvbnd;
                     rhsslack += intscalarrightrow * REALABS(valscurrentrow[c])
                        * REALABS(SCIPcolGetPrimsol(colscurrentrow[c]) - bestbndsol);
                     assert(SCIPisGE(scip, rhsslack, 0.0));
                     rhsisinfinity = rhsisinfinity || SCIPisInfinity(scip, rhs);                     
                  }  
               }

               if( !lhsisinfinity )
               {
                  if( valscurrentrow[c] * intscalarleftrow < 0 )
                  {
                     findClosestLb(scip, lpdata, colscurrentrow[c], &bestbndsol, &bestbndtype, &bestzvbnd, &bestbvbnd, &bestdvbnd );
                     assert(bestbndtype > -2 && lpdata->bestlbidxofcol[collppos] == bestbndtype);
                  }
                  else
                  {                  
                     assert(valscurrentrow[c] * intscalarleftrow > 0);
                     findClosestUb(scip, lpdata, colscurrentrow[c], &bestbndsol, &bestbndtype, &bestzvbnd, &bestbvbnd, &bestdvbnd );
                     assert(bestbndtype > -2 && lpdata->bestubidxofcol[collppos] == bestbndtype);
                  }
                  assert(bestbndtype == -1 || bestzvbnd != NULL);
                  
                  if( SCIPisInfinity(scip, -bestbndsol) || SCIPisInfinity(scip, bestbndsol) )
                     lhsisinfinity = TRUE;
                  else
                  {
                     /**@todo check whether REALABS is really correct */
                     if( bestbndtype == -1 )
                        lhs -= intscalarleftrow * REALABS(valscurrentrow[c]) * bestbndsol;
                     else
                        lhs -= intscalarleftrow * REALABS(valscurrentrow[c]) * bestdvbnd;
                     lhsslack += intscalarleftrow * REALABS(valscurrentrow[c])
                        * REALABS(SCIPcolGetPrimsol(colscurrentrow[c]) - bestbndsol);
                     assert(SCIPisGE(scip, lhsslack, 0.0));
                     lhsisinfinity = lhsisinfinity || SCIPisInfinity(scip, lhs);
                  }  
               }

               /* if both lhs and rhs became infinity, then the (relaxed) row is not relevant */
               if( lhsisinfinity && rhsisinfinity )
               {
                  rowisrelevant = FALSE;
                  break;
               }

               rowisrelevant = TRUE; 
            }

            /* check if row contains at least one relevant column (because k == 0 by initialization) */
            if( lpdata->subproblemsindexofcol[collppos] == k )
               rowisrelevant = TRUE;

            /* check if row contains no relevant columns but an odd lhs or rhs value */
            if( c == nnonzcurrentrow - 1 && (!lhsiseven || !rhsiseven) )
               rowisrelevant = TRUE;

         } 
         assert(SCIPisGE(scip, lhsslack, 0.0));
         assert(SCIPisGE(scip, rhsslack, 0.0));
         
         
         /* process row if it is relevant */
         if( rowisrelevant ) 
         {
            /* row is relevant because it contains a relevant column */
            problem->rrows[problem->nrrows] = r;
        
            lhsslackislessequalmaxslack = SCIPisLE(scip, lhsslack, maxslack);
            rhsslackislessequalmaxslack = SCIPisLE(scip, rhsslack, maxslack);

            /* note: due to the relaxation of continuous variables with their bounds the coeffs of nonzero variables
             * in left row and right row may be different. hence the row with smaller slack cannot be removed without
             * checking the coeffs first. 
             */
            if( !lhsisinfinity && lhsslackislessequalmaxslack )
            {
               /* "-a^T x <= -lhs" */
               lpdata->subproblemsindexofrow[r] = k;
               lpdata->rrowsindexofleftrow[r] = problem->nrrows;
              
               problem->rrows[problem->nrrows] = r;
               /**@todo check whether lhs is correct, or whether this must be -lhs. do we store the 
                * -ax <= -lhs constraint or the ax >= lhs constraint? Is this handled correctly above while updating lhs? 
                */ 
               problem->rrowsrhs[problem->nrrows] = lhs;
               problem->rrowsslack[problem->nrrows] = lhsslack;            
            }
            else
            {
               lpdata->rrowsindexofleftrow[r] = lhsisinfinity ? NONEXISTENT_ROW : SLACK_GREATER_THAN_MAXSLACK;
            }
            /**@todo check the following: if !lhsinfinity AND !rhsinfinity: then only one of them is stored currently, 
             * because problem->nrrows++ is only called once. if this is intended, why do we allocate 2 * lpdata->nrows 
             * entries for rrows?
             */
            if( !rhsisinfinity && rhsslackislessequalmaxslack )
            {
               /* "a^T x <= rhs" */
               lpdata->subproblemsindexofrow[r] = k;
               lpdata->rrowsindexofrightrow[r] = problem->nrrows;
              
               problem->rrows[problem->nrrows] = r;
               problem->rrowsrhs[problem->nrrows] = rhs;
               problem->rrowsslack[problem->nrrows] = rhsslack;    
            }
            else
            {
               lpdata->rrowsindexofrightrow[r] = rhsisinfinity ? NONEXISTENT_ROW : SLACK_GREATER_THAN_MAXSLACK;
            }

            /* increase counter only if at least one half row had a sufficiently small slack */
            if( lpdata->rrowsindexofleftrow[r] > -1 || lpdata->rrowsindexofrightrow[r] > -1 )            
               problem->nrrows++;
         }
         else /* case: !rowisrelevant */
         {
            /* row does not contain relevant columns */
            lpdata->subproblemsindexofrow[r] = IRRELEVANT;
            lpdata->rrowsindexofleftrow[r] = lhsisinfinity ? NONEXISTENT_ROW : NO_RELEVANT_COLUMNS;
            lpdata->rrowsindexofrightrow[r] = rhsisinfinity ? NONEXISTENT_ROW : NO_RELEVANT_COLUMNS;
         }      
      }
      else /* case: lhsslack > maxslack && rhsslack > maxslack */
      {
         lpdata->subproblemsindexofrow[r] = IRRELEVANT;
         lpdata->rrowsindexofleftrow[r] = lhsisinfinity ? NONEXISTENT_ROW : SLACK_GREATER_THAN_MAXSLACK;
         lpdata->rrowsindexofrightrow[r] = rhsisinfinity ? NONEXISTENT_ROW : SLACK_GREATER_THAN_MAXSLACK;
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &densecoeffscurrentleftrow);   
   SCIPfreeBufferArray(scip, &densecoeffscurrentrightrow);
   
   return SCIP_OKAY;
}


 
/* check if mod 2 data structure contains at most two nonzero entries per row */
static
SCIP_Bool hasMatrixMax2EntriesPerRow(
   ZEROHALF_MOD2DATA*    mod2data            /**< considered mod 2 data */
   )
{
   int r;
   int c;
   int nentries;

   assert(mod2data != NULL);

   if( mod2data->nrowsind == 0 )
      return TRUE; /*FALSE;*/

   if( mod2data->ncolsind <= 2 )
      return TRUE;

   for( r = 0; r < mod2data->nrowsind ; ++r )
   {
      nentries = 0;
      for( c = 0; c < mod2data->ncolsind ; ++c )
      {
         if( BITARRAYBITISSET(mod2data->rows[mod2data->rowsind[r]], mod2data->colsind[c]) ) /*lint !e701*/
         {
            nentries++;
            if( nentries > 2 )
               return FALSE;
         }
      }
   }
  
   return TRUE;
}


 
#ifdef ZEROHALF__PRINT_STATISTICS 
/* check if mod 2 data structure contains at most two nonzero entries per column */
static
SCIP_Bool hasMatrixMax2EntriesPerColumn(
   ZEROHALF_MOD2DATA*    mod2data            /**< considered mod 2 data */
   )
{
   int r;
   int c;
   int nentries;

   assert(mod2data != NULL);

   if( mod2data->ncolsind == 0 )
      return TRUE; /*FALSE;*/

   if( mod2data->nrowsind <= 2 )
      return TRUE;

   for( c = 0; c < mod2data->ncolsind ; ++c )
   {
      nentries = 0;
      for( r = 0; r < mod2data->nrowsind ; ++r )
      {
         if( BITARRAYBITISSET(mod2data->rows[mod2data->rowsind[r]], mod2data->colsind[c]) )
         {
            nentries++;
            if( nentries > 2 )
               return FALSE;
         }
      }
   }
  
   return TRUE;
}
#endif


 
/* stores relevant data into bit arrays (mod 2 data structure) */
static
SCIP_RETCODE storeMod2Data(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   int                   subproblemindex,    /**< index of considered subproblem */
   ZEROHALF_MOD2DATA*    mod2data            /**< data (mod 2) */
   )
{
   ZEROHALF_SUBLPDATA*   problem;
   SCIP_COL**            colscurrentrow;
   SCIP_ROW*             row;
   SCIP_Real*            nonzvalscurrentrow;
   SCIP_Real             maxslack;
   SCIP_Real             intscalar;
   BITARRAY              tempcurrentrow;   
   int*                  varboundstoadd;
   int                   nnonzcurrentrow;
   int                   rcolsindex;
   int                   c;
   int                   i;
   int                   j;
#ifdef ZEROHALF__PRINT_STATISTICS
   int                   nirrelevantvarbounds;
#endif
   SCIP_Bool             tempmod2rhs;
   SCIP_Bool             ignorerow;   
   SCIP_Bool             fliplhsrhs;
   SCIP_Bool             isrhsrow;
   SCIP_Real* densecoeffscurrentrow;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(lpdata != NULL);
   assert(lpdata->rows != NULL);
   assert(lpdata->nrows > 0);
   assert(lpdata->cols != NULL);
   assert(lpdata->ncols > 0);
   assert(lpdata->subproblems != NULL);
   assert(lpdata->nsubproblems > 0);
   assert(lpdata->subproblemsindexofrow != NULL);
   assert(lpdata->rrowsindexofleftrow != NULL);
   assert(lpdata->rrowsindexofrightrow != NULL);
   assert(lpdata->subproblemsindexofcol != NULL);
   assert(lpdata->rcolsindexofcol != NULL);
   assert(0 <= subproblemindex);
   assert(subproblemindex <= lpdata->nsubproblems);
   problem = lpdata->subproblems[subproblemindex];
   assert(problem != NULL);
   assert(problem->rrows != NULL);
   assert(problem->nrrows > 0);
   assert(problem->rrowsrhs != NULL);
   assert(problem->rrowsslack != NULL);
   assert(problem->rcols != NULL);
   assert(problem->nrcols > 0);
   assert(problem->rcolslbslack != NULL);
   assert(problem->rcolsubslack != NULL);
   assert(mod2data != NULL);
  
   /* identify varbounds to be added to the matrix */
   SCIP_CALL(SCIPallocBufferArray(scip, &varboundstoadd, 2 * problem->nrcols)); /* <0: lb, >0: ub */

   maxslack = sepadata->maxslack;

#ifdef ZEROHALF__PRINT_STATISTICS
   nirrelevantvarbounds = 0;
#endif
   mod2data->nvarbounds = 0;
   for( c = 0 ; c < problem->nrcols ; c++ )
   {
      SCIP_Bool lbslackisok;
      SCIP_Bool ubslackisok;

      lbslackisok = SCIPisLE(scip, problem->rcolslbslack[c], maxslack);
      ubslackisok = SCIPisLE(scip, problem->rcolsubslack[c], maxslack);

      if( lbslackisok && ubslackisok )
      {
         SCIP_Real lb;
         SCIP_Real ub;
         
         lb = SCIPcolGetLb(lpdata->cols[problem->rcols[c]]);
         ub = SCIPcolGetUb(lpdata->cols[problem->rcols[c]]);

         if( ISEVEN(scip, lb) != ISEVEN(scip, ub) )
         {
            varboundstoadd[mod2data->nvarbounds] = (-1) * (c + 1);
            mod2data->nvarbounds++;
            varboundstoadd[mod2data->nvarbounds] =  c + 1;
            mod2data->nvarbounds++;
         }
         else
         {
            if( SCIPisLE(scip, problem->rcolslbslack[c], problem->rcolsubslack[c]) )
               varboundstoadd[mod2data->nvarbounds] = (-1) * (c + 1);
            else
               varboundstoadd[mod2data->nvarbounds] =  c + 1;
            mod2data->nvarbounds++;
#ifdef ZEROHALF__PRINT_STATISTICS
            nirrelevantvarbounds++;
#endif
         }
      }
      else
      {
         if( lbslackisok )
         {
            varboundstoadd[mod2data->nvarbounds] = (-1) * (c + 1);
            mod2data->nvarbounds++;
         }
#ifdef ZEROHALF__PRINT_STATISTICS
         else
            nirrelevantvarbounds++;
#endif
         if( ubslackisok )
         {
            varboundstoadd[mod2data->nvarbounds] =  c + 1;
            mod2data->nvarbounds++;
         }
#ifdef ZEROHALF__PRINT_STATISTICS
         else
            nirrelevantvarbounds++;
#endif
      }
   }
   mod2data->nrows = problem->nrrows + mod2data->nvarbounds;
  
   /* allocate temporary memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mod2data->rows), mod2data->nrows) ); 
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mod2data->rowaggregations), mod2data->nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mod2data->rhs), mod2data->nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mod2data->slacks), mod2data->nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mod2data->fracsol), problem->nrcols) ); 
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mod2data->rowstatistics), mod2data->nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mod2data->colstatistics), problem->nrcols) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mod2data->rowsind), mod2data->nrows) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mod2data->colsind), problem->nrcols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &densecoeffscurrentrow, lpdata->ncols) );
  
   /* initialize temporary memory */
   mod2data->relatedsubproblem = problem;
   BMSclearMemoryArray(mod2data->rows, mod2data->nrows);                 /* NULL = 0x0 */
   BMSclearMemoryArray(mod2data->rowaggregations, mod2data->nrows);      /* NULL = 0x0 */
   BMScopyMemoryArray(mod2data->slacks, problem->rrowsslack, problem->nrrows); 
   for( c = 0 ; c < problem->nrcols ; ++c)
      mod2data->fracsol[c] = SCIPcolGetPrimsol(lpdata->cols[problem->rcols[c]]);
   for( c = 0 ; c < problem->nrcols ; c++)
      mod2data->colsind[c] = c;
   mod2data->nrowsind = 0;
   mod2data->ncolsind = problem->nrcols;
   mod2data->rowsbitarraysize = (int) GETREQUIREDBITARRAYSIZE(problem->nrcols); 
   mod2data->rowaggregationsbitarraysize = (int) GETREQUIREDBITARRAYSIZE(problem->nrrows);
   tempcurrentrow = NULL;

   /* (i) for all relevant rows */
   for( i = 0 ; i < problem->nrrows ; ++i )
   {
      row = lpdata->rows[problem->rrows[i]]; 
      colscurrentrow = SCIProwGetCols(row);
      nonzvalscurrentrow = SCIProwGetVals(row);
      nnonzcurrentrow = SCIProwGetNLPNonz(row);
      assert(nnonzcurrentrow > 0);
      tempcurrentrow = NULL;
      fliplhsrhs = FALSE;
      ignorerow = FALSE;      

      /* check if rrows corresponds to a lhs or rhs row in the LP */           
      if( lpdata->rrowsindexofleftrow[problem->rrows[i]] == i )
         isrhsrow = FALSE;
      else
      {
         assert(lpdata->rrowsindexofrightrow[problem->rrows[i]] == i);
         isrhsrow = TRUE;             
      }
      intscalar = isrhsrow ? lpdata->intscalarsrightrow[problem->rrows[i]]
         : lpdata->intscalarsleftrow[problem->rrows[i]]; 
    
      /* clear dense coeffs array */
      BMSclearMemoryArray(densecoeffscurrentrow, lpdata->ncols);

      /* compute dense coeffs array of current row (including intscaling and bound substitutions) */
      for( j = 0 ; j < nnonzcurrentrow; ++j )
      {
         if( SCIPvarGetType(SCIPcolGetVar(colscurrentrow[j])) == SCIP_VARTYPE_CONTINUOUS )
         {    
            SCIP_Bool ispositivecoeff;           
            SCIP_VAR* bestzvbnd;
            SCIP_Real bestbndsol;
            SCIP_Real bestbvbnd;
            SCIP_Real bestdvbnd;
            int bestbndtype;

            /* check sign of coefficient */
            if( nonzvalscurrentrow[j] * intscalar > 0.0 )
               ispositivecoeff = TRUE;
            else
               ispositivecoeff = FALSE;
            
            /* get appropriate bound */
            if( isrhsrow == ispositivecoeff )           
               findClosestLb(scip, lpdata, colscurrentrow[j], &bestbndsol, &bestbndtype, &bestzvbnd, &bestbvbnd, &bestdvbnd );
            else
               findClosestUb(scip, lpdata, colscurrentrow[j], &bestbndsol, &bestbndtype, &bestzvbnd, &bestbvbnd, &bestdvbnd );

            /* check bound type */
            assert(bestbndtype > -2);
            if( !USEVARBOUNDS )  /*lint !e774 !e506*/
               assert(bestbndtype == -1);

            /* normal lb or ub is used; only rhs would have to be adjusted but this has already been done in getRelevantRows */
            if( bestbndtype == -1 )
               continue;
            assert(USEVARBOUNDS && bestbndtype > -1); /*lint !e774 !e506*/

            /* variable bound is used: update coefficient of non-continuous variable z that is used in substitution */
            densecoeffscurrentrow[SCIPcolGetLPPos(SCIPvarGetCol(bestzvbnd))] += (nonzvalscurrentrow[j] * intscalar * bestbvbnd);
         }         
         else
         {
            densecoeffscurrentrow[SCIPcolGetLPPos(colscurrentrow[j])] += (nonzvalscurrentrow[j] * intscalar);
         }
      }         
      
      for( j = 0 ; j < lpdata->ncols; ++j )
      {
         assert(SCIPcolGetLPPos(lpdata->cols[j]) == j);

         if( SCIPisZero(scip, densecoeffscurrentrow[j]) )
            continue;
         
         if( intscalar == 1.0 && !SCIPisIntegral(scip, densecoeffscurrentrow[j]) )
         {
            ignorerow = TRUE;
            break;
         }
         else
            assert(sepadata->scalefraccoeffs);

         /* integral coefficient */
         /* coefficient is only integral with respect to tolerances; use really integral values */
         if( isrhsrow )
            densecoeffscurrentrow[j] = SCIPfloor(scip, densecoeffscurrentrow[j]);
         else
            densecoeffscurrentrow[j] = SCIPceil(scip, densecoeffscurrentrow[j]);

         if( ISODD(scip, densecoeffscurrentrow[j]) ) 
         {
            rcolsindex = lpdata->rcolsindexofcol[j];
            fliplhsrhs = XOR((int) fliplhsrhs, 
               (int) (rcolsindex == LP_SOL_EQUALS_ODD_LB || rcolsindex == LP_SOL_EQUALS_ODD_UB));
            if( rcolsindex >= 0 ) /* relevant column? */
            {
               if( tempcurrentrow == NULL )
               {
                  SCIP_CALL( SCIPallocMemoryArray(scip, &tempcurrentrow, mod2data->rowsbitarraysize) );
                  BITARRAYCLEAR(tempcurrentrow, mod2data->rowsbitarraysize);
               }
               assert(rcolsindex < problem->nrcols);
               BITARRAYBITSET(tempcurrentrow, rcolsindex); /*lint !e701*/
               assert(BITARRAYBITISSET(tempcurrentrow, rcolsindex)); /*lint !e701*/
            }
         }
      }

      /* check if current row should be ignored, continuing with the next one */
      if( ignorerow )
      {
         if( tempcurrentrow != NULL )
         {
            SCIPfreeMemoryArray(scip, &tempcurrentrow);
            tempcurrentrow = NULL;
         }
         continue;
      }
    
      /* consider rhs */
      if( XOR((int) ISODD(scip, problem->rrowsrhs[i]), (int) fliplhsrhs) )
         tempmod2rhs = TRUE;
      else
         tempmod2rhs = FALSE;
    
      if( tempcurrentrow == NULL && tempmod2rhs )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &tempcurrentrow, mod2data->rowsbitarraysize) );
         BITARRAYCLEAR(tempcurrentrow, mod2data->rowsbitarraysize);
      }
      assert(tempcurrentrow != NULL || !tempmod2rhs);
    
      /* store temporary data in appropriate (mod 2) data structures */
      if( tempcurrentrow != NULL )
      { 
         mod2data->rows[i] = tempcurrentrow;
         mod2data->rhs[i] = tempmod2rhs;

         assert(mod2data->rowaggregationsbitarraysize > 0);
         SCIP_CALL( SCIPallocMemoryArray(scip, &(mod2data->rowaggregations[i]), 
               mod2data->rowaggregationsbitarraysize) ); /*lint !e866*/
         BITARRAYCLEAR(mod2data->rowaggregations[i], mod2data->rowaggregationsbitarraysize); /*lint !e866*/
         BITARRAYBITSET(mod2data->rowaggregations[i], i); /*lint !e701*/

         mod2data->rowsind[mod2data->nrowsind] = i;
         mod2data->nrowsind++;
      
         tempcurrentrow = NULL;
      }
      else
      {
         /* zero row */
         lpdata->subproblemsindexofrow[problem->rrows[i]] = IRRELEVANT;
         if( lpdata->rrowsindexofleftrow[problem->rrows[i]] >= 0 )
            lpdata->rrowsindexofleftrow[problem->rrows[i]] = ZERO_ROW;
         if( lpdata->rrowsindexofrightrow[problem->rrows[i]] >= 0 )
            lpdata->rrowsindexofrightrow[problem->rrows[i]] = ZERO_ROW;
      }  
   }

  
   /* (ii)   for all relevant varbounds */
   i = problem->nrrows;
   for( j = 0 ; j < mod2data->nvarbounds ; ++j)
   {
      SCIP_Real bound;

      if( varboundstoadd[j] < 0 )
         c = (-1) * varboundstoadd[j] - 1;
      else
         c = varboundstoadd[j] - 1;

      assert(mod2data->rowsbitarraysize > 0);
      SCIP_CALL(SCIPallocMemoryArray(scip, &(mod2data->rows[i]), mod2data->rowsbitarraysize)); /*lint !e866*/
      BITARRAYCLEAR(mod2data->rows[i], mod2data->rowsbitarraysize); /*lint !e866*/
      BITARRAYBITSET(mod2data->rows[i], c); /*lint !e701*/
      assert(BITARRAYBITISSET(mod2data->rows[i], c)); /*lint !e701*/

      SCIP_CALL(SCIPallocMemoryArray(scip, &(mod2data->rowaggregations[i]), 
            mod2data->rowaggregationsbitarraysize)); /*lint !e866*/
      BITARRAYCLEAR(mod2data->rowaggregations[i], mod2data->rowaggregationsbitarraysize); /*lint !e866*/
    
      if( varboundstoadd[j] < 0 )
      {
         bound = SCIPcolGetLb(lpdata->cols[problem->rcols[c]]);
         mod2data->rhs[i] = ISODD(scip, bound);
         mod2data->slacks[i] = problem->rcolslbslack[c];
      }
      else
      {
         bound = SCIPcolGetUb(lpdata->cols[problem->rcols[c]]);
         mod2data->rhs[i] = ISODD(scip, bound);
         mod2data->slacks[i] = problem->rcolsubslack[c];
      }
      if( SCIPisFeasZero(scip, mod2data->slacks[i]) )
         mod2data->slacks[i] = 0.0;

      mod2data->rowsind[mod2data->nrowsind] = i;
      mod2data->nrowsind++;
      i++;
   }
  
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &densecoeffscurrentrow);   
   SCIPfreeBufferArray(scip, &varboundstoadd); 
  
#ifdef ZEROHALF__PRINT_STATISTICS
   ZEROHALFstatisticsMessage("\n");
   ZEROHALFstatisticsMessage("                | ------------------------------- subproblem ------------------------------- | ------------------------------\n");
   ZEROHALFstatisticsMessage("                |   nrrows |   nrcols | nvarbnds | ndlvbnds | max2/row | max2/col |  A^T ept |                               \n");
   ZEROHALFstatisticsMessage("%15s | %8d | %8d | %8d | %8d | %8s | %8s | %8s |\n",
      "SUBPROBLEMDATA", problem->nrrows, problem->nrcols, mod2data->nvarbounds, nirrelevantvarbounds,
      hasMatrixMax2EntriesPerRow(mod2data) ? "yes" : "no", hasMatrixMax2EntriesPerColumn(mod2data) ? "yes" : "no", "n/a");
#endif
  
   return SCIP_OKAY;
}


/* --------------------------------------------------------------------------------------------------------------------
 * local methods: cut generation
 * -------------------------------------------------------------------------------------------------------------------- */


/** stores nonzero elements of dense coefficient vector as sparse vector, and calculates activity and norm */
static
SCIP_RETCODE storeCutInArrays(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of problem variables */
   SCIP_VAR**            vars,               /**< problem variables */
   SCIP_Real*            cutcoefs,           /**< dense coefficient vector */
   SCIP_Real*            varsolvals,         /**< dense variable LP solution vector */
   char                  normtype,           /**< type of norm to use for efficacy norm calculation */

   SCIP_VAR**            cutvars,            /**< array to store variables of sparse cut vector */
   SCIP_Real*            cutvals,            /**< array to store coefficients of sparse cut vector */
   int*                  cutlen,             /**< pointer to store number of nonzero entries in cut */
   SCIP_Real*            cutact,             /**< pointer to store activity of cut */
   SCIP_Real*            cutnorm             /**< pointer to store norm of cut vector */
   )
{
   SCIP_Real             val;
   SCIP_Real             absval;
   SCIP_Real             cutsqrnorm;
   SCIP_Real             act;
   SCIP_Real             norm;
   int                   len;
   int                   v;

   assert(nvars == 0 || cutcoefs != NULL);
   assert(nvars == 0 || varsolvals != NULL);
   assert(cutvars != NULL);
   assert(cutvals != NULL);
   assert(cutlen != NULL);
   assert(cutact != NULL);
   assert(cutnorm != NULL);

   len = 0;
   act = 0.0;
   norm = 0.0;
   switch(normtype)
   {
   case 'e':
      cutsqrnorm = 0.0;
      for( v = 0; v < nvars; ++v)
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            cutsqrnorm += SQR(val);
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      norm = SQRT(cutsqrnorm);
      break;
   case 'm':
      for( v = 0; v < nvars; ++v)
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            absval = REALABS(val);
            norm = MAX(norm, absval);
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      break;
   case 's':
      for( v = 0; v < nvars; ++v)
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            norm += REALABS(val);
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      break;
   case 'd':
      for( v = 0; v < nvars; ++v)
      {
         val = cutcoefs[v];
         if( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            norm = 1.0;
            cutvars[len] = vars[v];
            cutvals[len] = val;
            len++;
         }
      }
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", normtype);
      return SCIP_INVALIDDATA;
   }

   *cutlen = len;
   *cutact = act;
   *cutnorm = norm;

   return SCIP_OKAY;
}


 
/** adds a separated zerohalf cut to SCIP if it was successfully created and is efficacious */
static
SCIP_RETCODE addZerohalfCutToLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   ZEROHALF_CUTDATA*     cutdata,            /**< separated zerohalf cut */
   int*                  nsepacuts,          /**< pointer to store number of separated (efficacious) zerohalf cuts */
   SCIP_RESULT*          result              /**< pointer to store return code */
   )
{
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(cutdata != NULL);
   assert(result != NULL);

   /* check if SCIPcalcMIR was not successful */
   if( !cutdata->isfeasviolated || !cutdata->success )
      return SCIP_OKAY;

   /* check if norm was not calculated correctly */
   if( !SCIPisPositive(scip, cutdata->norm) )
   {
      SCIPerrorMessage("Zerohalf cut norm is NOT positive!\n");    
      return SCIP_ERROR;
   }

   /* check if cut is not efficacious */
   if( !sepadata->forcecutstolp && !sepadata->forcecutstosepastore
      && !SCIPisEfficacious(scip, cutdata->efficacy) )
   {
      return SCIP_OKAY;
   }
          
   /* add cut (if no cutpool is used otherwise add it at the end of the separation main method)*/
   if( !sepadata->usezhcutpool )
   {
      SCIP_CALL(SCIPaddCut(scip, NULL, cutdata->cut, sepadata->forcecutstolp));
      if( !cutdata->islocal )
      {
         SCIP_CALL(SCIPaddPoolCut(scip, cutdata->cut));
      }
   }
  
   cutdata->addedtolp = TRUE;
   (*nsepacuts)++;           

   *result = SCIP_SEPARATED;
          
   SCIPdebug( SCIP_CALL( SCIPprintRow(scip, cutdata->cut, NULL) ) );

   return SCIP_OKAY;
}


/* --------------------------------------------------------------------------------------------------------------------
 * local methods: preprocessing
 * -------------------------------------------------------------------------------------------------------------------- */


 
/** marks a row as "removed" and stores why it has been removed using a flag */
static
void markRowAsRemoved(
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered mod 2 data */
   int                   r,                  /**< mod2data->rows index of row that shall be removed */
   int                   flag                /**< flag (cause of removal) */
   )
{
   assert(mod2data != NULL);
   assert(mod2data->relatedsubproblem != NULL);
   assert(mod2data->nrowsind > 0);
   assert(r >= 0);
   assert(r < mod2data->nrowsind);

   mod2data->rowstatistics[mod2data->rowsind[r]] = flag;
}


 
/** marks a row as "removed" and stores why it has been removed using a flag. in addition it clears this column's mod 2 data */
static
void  markColAsRemovedAndClearCol(
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered mod 2 data */
   int                   c,                  /**< mod2data->relatedsubproblem->rcols index of column that shall be removed */
   int                   flag                /**< flag (cause of removal) */
   )
{
   int                   i;
   int                   rowsbind;
   BITARRAYBITMASKTYPE   rowsbmask;

   assert(mod2data != NULL);
   assert(mod2data->relatedsubproblem != NULL);
   assert(mod2data->ncolsind > 0);
   assert(c >= 0);
   assert(c < mod2data->ncolsind);

  
   /* mark col */
   mod2data->colstatistics[mod2data->colsind[c]] = flag;

   /* clear col */
   rowsbind = (int) GETBITARRAYINDEX(mod2data->colsind[c]);
   rowsbmask = ~GETBITARRAYMASK(mod2data->colsind[c]); /*lint !e701*/
   for( i = 0 ; i < mod2data->nrowsind ; ++i)
      mod2data->rows[mod2data->rowsind[i]][rowsbind] &= rowsbmask;  
}

 


/** given a subset of mod 2 rows it returns a {0,1/2} weight vector used to
    combine the (original) LP rows. Note: original rows a stored as lhs <= a^Tx
    <= rhs by SCIP. Positive weights refer to "right half-rows" a^Tx <= rhs and
    negative weights to "left half-rows" -a^Tx <= -lhs */ 
static
SCIP_RETCODE getZerohalfWeightvectorFromSelectedRowsBitarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */ 
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered mod 2 data */
   BITARRAY              rrowsincut,         /**< subset of selected mod2data->rows */
   SCIP_Real**           weights,            /**< pointer to store the {-0.5,0,0.5} weights vector */
   int*                  nrowsincut          /**< pointer to store the number of combined original LP rows */   
   )
{   /*lint --e{438}*/

   ZEROHALF_SUBLPDATA*   problem;
   int                   lppos;
   int                   i;
   int                   nnonz;

   assert(scip != NULL);
   assert(lpdata != NULL);
   assert(lpdata->nrows > 0);
   assert(mod2data != NULL);
   assert(mod2data->relatedsubproblem != NULL);
   assert(mod2data->nrowsind > 0);
   assert(rrowsincut != NULL);
   assert(weights != NULL);
   assert(*weights == NULL);
   assert(nrowsincut != NULL);

  
   /* allocate temporary memory */ 
   SCIP_CALL(SCIPallocMemoryArray(scip, weights, lpdata->nrows));

   /* initialize */
   BMSclearMemoryArray(*weights, lpdata->nrows);
   problem = mod2data->relatedsubproblem;
  
   /* determine row weights */
   *nrowsincut = 0;
   nnonz = 0;
   for( i = 0 ; i < problem->nrrows ; ++i)
   {
      lppos = problem->rrows[i];
      assert(0 <= lppos && lppos <= lpdata->nrows);         
      if( BITARRAYBITISSET(rrowsincut, i) ) /*lint !e701*/
      {
         assert(lpdata->rrowsindexofleftrow[lppos] == i || lpdata->rrowsindexofrightrow[lppos] == i);
         
         SCIPdebugMessage("  %1s0.5   (int scaling: %16.4f / %16.4f)  row[%d] %s\n",
            lpdata->rrowsindexofleftrow[lppos] == i ? "-" : "+",
            lpdata->intscalarsleftrow[lppos], lpdata->intscalarsrightrow[lppos],
            lppos, SCIProwGetName(lpdata->rows[lppos]));          
      
         if( lpdata->rrowsindexofleftrow[lppos] == i )
            (*weights)[lppos] = lpdata->intscalarsleftrow[lppos] * (-0.5);
         else
            (*weights)[lppos] = lpdata->intscalarsrightrow[lppos] * 0.5;
      
         nnonz += SCIProwGetNLPNonz(lpdata->rows[lppos]); 
         (*nrowsincut)++;
      }
   }

   /* check if row aggregation might be too dense */
   if( nnonz >= 5 * sepadata->maxnnonz ) 
   {
      SCIPfreeMemoryArray(scip, weights);
      weights = NULL;/*lint !e438*/
   }

   return SCIP_OKAY;
}


/** creates a zerohalf cut from a given weightvector */
static
SCIP_RETCODE createZerohalfCutFromZerohalfWeightvector(
   SCIP*                 scip,               /**< SCIP data structure */              
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   SCIP_Real*            weights,            /**< weightvector */
   char                  normtype,           /**< SCIP normtype */
   int                   nzerohalfcuts,      /**< number of zerohalf cuts (used for naming the cut) */
   SCIP_Real**           varsolvals,         /**< pointer to array of LP solution values of variables */
   ZEROHALF_CUTDATA*     cutdata             /**< pointer to data structure used for storing the cut */
   )
{
   SCIP_Real*            cutcoefs;
   SCIP_VAR**            cutvars;
   SCIP_Real*            cutvals;
   char                  cutname[SCIP_MAXSTRLEN];
  
   assert(scip != NULL);
   assert(lpdata != NULL);
   assert(lpdata->nvars > 0);
   assert(weights != NULL);
   assert(varsolvals != NULL); 
   assert(cutdata != NULL);
   assert(cutdata->relatedsubproblem != NULL);

   /* note: cutdata->relatedmod2data can be NULL if cut was determined
    *       before mod 2 data structures were created */
  
   /* allocate temporary memory */
   SCIP_CALL(SCIPallocBufferArray(scip, &cutcoefs, lpdata->nvars));
  
   /* calculate MIR */
   cutdata->success = FALSE;
   if( sepadata->maxtestdelta == 0 )
   {
      /* generate cut for delta = 1.0 */
      SCIP_CALL( SCIPcalcMIR(scip, NULL, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, FIXINTEGRALRHS,
            BOUNDSFORTRANS, BOUNDTYPESFORTRANS, sepadata->maxnnonz, MAXWEIGHTRANGE, MINFRAC, MAXFRAC,
            weights, 1.0, NULL, NULL, cutcoefs, &(cutdata->rhs), &(cutdata->activity),
            &(cutdata->success), &(cutdata->islocal)) );
     
      if( sepadata->trynegscaling )
      {
         SCIP_CALL( SCIPcalcMIR(scip, NULL, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, FIXINTEGRALRHS,
               BOUNDSFORTRANS, BOUNDTYPESFORTRANS, sepadata->maxnnonz, MAXWEIGHTRANGE, MINFRAC, MAXFRAC,
               weights, -1.0, NULL, NULL, cutcoefs, &(cutdata->rhs), &(cutdata->activity),
               &(cutdata->success), &(cutdata->islocal)) );
      }
   }
   else
   {
      int ncuts; 
      SCIP_Real bestdelta;
      SCIP_Bool bestdeltavalid;

      ncuts = 0; 

      if( *varsolvals == NULL )
      {
         /* get the solution values for all active variables */
         SCIP_CALL(SCIPallocMemoryArray(scip, varsolvals, lpdata->nvars));
         SCIP_CALL( SCIPgetSolVals(scip, NULL, lpdata->nvars, lpdata->vars, *varsolvals) );

#ifndef NDEBUG
         /* because later when calling SCIPcutGenerationHeuristicCmir() varsolvals are used, it is needed that the
          * corresponding variables have the same order here and there, so we do the same checking and test that all
          * variables are ordered by their problem index
          */
         {
            int i;
            for(i = lpdata->nvars - 1; i >= 0; --i )
               assert(i == SCIPvarGetProbindex(lpdata->vars[i]));
         }
#endif
      }
      assert(*varsolvals != NULL);

      /* find best value of delta */
      SCIP_CALL( SCIPcutGenerationHeuristicCmir(scip, sepa, NULL, *varsolvals, sepadata->maxtestdelta, weights, BOUNDSWITCH, 
            USEVBDS, ALLOWLOCAL, FIXINTEGRALRHS, sepadata->maxnnonz, MAXWEIGHTRANGE, MINFRAC, MAXFRAC, 
            sepadata->trynegscaling, TRUE, "zerohalf", &ncuts, &bestdelta, &bestdeltavalid) );  
      assert(ncuts == 0);

      /* best delta corresponds to an efficient cut */
      if( bestdeltavalid ) 
      {  
         SCIP_CALL( SCIPcalcMIR(scip, NULL, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, FIXINTEGRALRHS,
               BOUNDSFORTRANS, BOUNDTYPESFORTRANS, sepadata->maxnnonz, MAXWEIGHTRANGE, MINFRAC, MAXFRAC,
               weights, bestdelta, NULL, NULL, cutcoefs, &(cutdata->rhs), &(cutdata->activity),
               &(cutdata->success), &(cutdata->islocal)) );
      }
   }
   assert(ALLOWLOCAL || !cutdata->islocal);
  
   cutdata->violation = cutdata->activity - cutdata->rhs;
  
   /* if successful, convert dense cut into sparse row */
   if( cutdata->success )
   {
      cutdata->isfeasviolated = SCIPisFeasGT(scip, cutdata->activity, cutdata->rhs);
      SCIPdebugMessage("Cut is %sfeasviolated: (act: %e, rhs: %e, viol: %e)\n", 
         cutdata->isfeasviolated ? "" : "not ", cutdata->activity, cutdata->rhs, cutdata->violation);      

      if( cutdata->isfeasviolated )    
      { 
         if( *varsolvals == NULL )
         {
            /* get the solution values for all active variables */
            SCIP_CALL(SCIPallocMemoryArray(scip, varsolvals, lpdata->nvars));
            SCIP_CALL( SCIPgetSolVals(scip, NULL, lpdata->nvars, lpdata->vars, *varsolvals) );
         }
         assert(*varsolvals != NULL);
      
         /* get temporary memory for storing the cut as sparse row */      
         SCIP_CALL(SCIPallocBufferArray(scip, &cutvars, lpdata->nvars));
         SCIP_CALL(SCIPallocBufferArray(scip, &cutvals, lpdata->nvars));
      
         /* store the cut as sparse row, calculate activity and norm of cut */
         SCIP_CALL(storeCutInArrays(scip, lpdata->nvars, lpdata->vars,
               cutcoefs, *varsolvals, normtype, cutvars, cutvals,
               &(cutdata->nnonz), &(cutdata->activity), &(cutdata->norm)));

         /* check cut norm and efficacy */
         if( SCIPisPositive(scip, cutdata->norm) )
         {
            cutdata->efficacy = (cutdata->activity - cutdata->rhs) / cutdata->norm;

            if( sepadata->forcecutstolp || sepadata->forcecutstosepastore
               || (SCIPisEfficacious(scip, cutdata->efficacy) && cutdata->nnonz < sepadata->maxnnonz) )
            {    
               /* create cut */
               (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN,"zerohalf%d_%d", SCIPgetNLPs(scip), nzerohalfcuts);
               SCIP_CALL(SCIPcreateEmptyRowSepa(scip, &(cutdata->cut), sepa, cutname, -SCIPinfinity(scip), cutdata->rhs, 
                     cutdata->islocal, FALSE, sepadata->dynamiccuts));
               SCIP_CALL(SCIPaddVarsToRow(scip, cutdata->cut, cutdata->nnonz, cutvars, cutvals));
            }
            else
               cutdata->success = FALSE;
         }
         else
            cutdata->success = FALSE;

         /* free temporary memory */
         SCIPfreeBufferArray(scip, &cutvals);
         SCIPfreeBufferArray(scip, &cutvars);
      }
      else
         cutdata->success = FALSE;
   } 

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &cutcoefs);
  
   return SCIP_OKAY;
}


/** searches for trivial zerohalf cuts, given as (0,..0) row with rhs=1 and slack <= maxslack */
static
SCIP_RETCODE preprocessTrivialZerohalfCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */        
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered (preprocessed) subproblem mod 2 */
   int                   firstrowsind,       /**< first mod2data->rows index to be considered */
   int                   lastrowsind,        /**< last mod2data->rows index to be considered */
   char                  normtype,           /**< SCIP normtype */
   int                   maxsepacuts,        /**< maximal number of zerohalf cuts separated per separation round */
   int                   maxcuts,            /**< maximal number of zerohalf cuts found per separation round (incl. ineff. cuts) */
   int*                  nsepacuts,          /**< pointer to store current number of separated zerohalf cuts */
   int*                  nzerohalfcuts,      /**< pointer to store current number of found zerohalf cuts */
   ZEROHALF_CUTDATA**    zerohalfcuts,       /**< pointer to store a found zerohalf cut */
   SCIP_Real**           varsolvals,         /**< dense variable LP solution vector */
   CUTSEPARATEDBY        cutseparatedby,     /**< flag */
   SCIP_RESULT*          result              /**< pointer to SCIP result value of separation */
   )
{
   int                   r;
   int                   r2;
   int                   nrowsremoved;
   SCIP_Real             maxslack;
   BITARRAY              zerorow;
   SCIP_Bool*            removerow;
   SCIP_Real*            weights;
   int                   nrowsincut;
    
   assert(scip != NULL);
   assert(lpdata != NULL);
   assert(mod2data != NULL);
   assert(mod2data->relatedsubproblem != NULL);
   assert(firstrowsind >= 0);
   assert(lastrowsind <= mod2data->nrowsind);
   assert(nsepacuts != NULL);
   assert(nzerohalfcuts != NULL);
   assert(zerohalfcuts != NULL);
   assert(*nsepacuts <= *nzerohalfcuts);
   assert(varsolvals != NULL);
   assert(result != NULL);
  
  
   /* check if matrix or colind range is empty */
   if( mod2data->nrowsind == 0 || lastrowsind - firstrowsind <= 0 )
      return SCIP_OKAY;


   /* allocate temporary memory */
   SCIP_CALL(SCIPallocBufferArray(scip, &removerow, lastrowsind - firstrowsind));
   SCIP_CALL(SCIPallocBufferArray(scip, &zerorow, mod2data->rowsbitarraysize));
  
   /* initialize */
   BMSclearMemoryArray(zerorow, mod2data->rowsbitarraysize);
   BMSclearMemoryArray(removerow, lastrowsind - firstrowsind);  
   maxslack = sepadata->maxslack;
   nrowsremoved = 0;

  
   /* check all rows */
   for( r = 0 ; r < lastrowsind - firstrowsind && *nsepacuts < maxsepacuts && *nzerohalfcuts < maxcuts; ++r )
      if( mod2data->rhs[mod2data->rowsind[firstrowsind + r]] == TRUE )
         if( SCIPisLE(scip, mod2data->slacks[mod2data->rowsind[firstrowsind + r]], maxslack ))
            if( BITARRAYSAREEQUAL(mod2data->rows[mod2data->rowsind[firstrowsind + r]],
                  zerorow, mod2data->rowsbitarraysize) ) /* check if row is (0 ... 0 , 1) */
            {
               /* a violated zerohalf cut has been found */
               weights = NULL;
               SCIP_CALL(getZerohalfWeightvectorFromSelectedRowsBitarray(scip, sepadata, lpdata, mod2data,
                     mod2data->rowaggregations[mod2data->rowsind[firstrowsind + r]], &weights, &nrowsincut));
               if( weights == NULL )
               {
                  continue;
               }
               assert(nrowsincut > 0);

               /* create zerohalf cut */
               SCIP_CALL(ZerohalfCutDataCreate(scip, &(zerohalfcuts[*nzerohalfcuts]),
                     mod2data->relatedsubproblem, mod2data, 1, nrowsincut, cutseparatedby));
               SCIP_CALL(createZerohalfCutFromZerohalfWeightvector(scip, sepa, sepadata,
                     lpdata, weights, normtype, *nzerohalfcuts, varsolvals, zerohalfcuts[*nzerohalfcuts]));
                    
               /* add cut */
               SCIP_CALL( addZerohalfCutToLP(scip, sepadata, zerohalfcuts[*nzerohalfcuts], nsepacuts, result) );
               (*nzerohalfcuts)++;
          
               /* free temporary memory */
               SCIPfreeMemoryArray(scip, &weights);

               removerow[r] = TRUE;
               nrowsremoved++;
            }


   /* update mod2data->rowsind if necessary */
   if( nrowsremoved > 0 )
   {   
      r2 = firstrowsind;
      for( r = firstrowsind ; r < mod2data->nrowsind && r2 < mod2data->nrowsind; ++r)
      {
         if( r < lastrowsind - firstrowsind )
            while( removerow[r2] && r2 < mod2data->nrowsind )
               r2++;
         if( r < r2 && r2 < mod2data->nrowsind )
            mod2data->rowsind[r] = mod2data->rowsind[r2];
         r2++;
      }      
      mod2data->nrowsind -= nrowsremoved;
   }
  
  
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &zerorow);
   SCIPfreeBufferArray(scip, &removerow);

   return SCIP_OKAY;
}


/** applies some row reductions */
static
SCIP_RETCODE preprocessRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */       
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered (preprocessed) subproblem mod 2 */
   int                   firstrowsind,       /**< first mod2data->rows index to be considered */
   int                   lastrowsind,        /**< last mod2data->rows index to be considered */
   SCIP_Bool             removezerorows,     /**< should zero rows be removed? */
   SCIP_Bool             removelargeslackrows, /**< should rows with slack > maxslack be removed? */
   SCIP_Bool             removeidenticalrows /**< should identical rows be removed? */
   )
{
   int                   r1;
   int                   r2;
   SCIP_Bool*            rowisprocessed;
   SCIP_Bool*            removerow;
   int                   nzerorowsremoved;
   int                   nlargeslackrowsremoved;
   int                   nidenticalrowsremoved;
   SCIP_Real             maxslack;
   BITARRAY              zerorow;
  
   assert(scip != NULL);
   assert(lpdata != NULL);
   assert(mod2data != NULL);
   assert(mod2data->relatedsubproblem != NULL);
   assert(firstrowsind >= 0);
   assert(lastrowsind <= mod2data->nrowsind);

  
   /* check if matrix or colind range is empty */
   if( mod2data->nrowsind == 0 || lastrowsind - firstrowsind <= 0 )
      return SCIP_OKAY;

   /* allocate temporary memory */
   SCIP_CALL(SCIPallocBufferArray(scip, &rowisprocessed, lastrowsind - firstrowsind));
   SCIP_CALL(SCIPallocBufferArray(scip, &removerow, lastrowsind - firstrowsind));
   zerorow = NULL;
   if( removezerorows )
   {
      SCIP_CALL(SCIPallocBufferArray(scip, &zerorow, mod2data->rowsbitarraysize));
   }
   /* initialize */
   BMSclearMemoryArray(rowisprocessed, lastrowsind - firstrowsind);
   BMSclearMemoryArray(removerow, lastrowsind - firstrowsind);  
   if( removezerorows )
   {
      BMSclearMemoryArray(zerorow, mod2data->rowsbitarraysize);    
   }
   maxslack = sepadata->maxslack;
   nzerorowsremoved = 0;
   nlargeslackrowsremoved = 0;
   nidenticalrowsremoved = 0;

   /* check all pairs of rows */
   for( r1 = 0 ; r1 < lastrowsind - firstrowsind ; ++r1)
      if( !rowisprocessed[r1] )
      {
         rowisprocessed[r1] = TRUE;

         if( removezerorows && !removerow[r1] )
            if( mod2data->rhs[mod2data->rowsind[firstrowsind + r1]] == FALSE )
            {
               assert(zerorow != NULL);
               if( BITARRAYSAREEQUAL(mod2data->rows[mod2data->rowsind[firstrowsind + r1]],
                     zerorow, mod2data->rowsbitarraysize) )
               {
                  markRowAsRemoved(mod2data, firstrowsind + r1, ZERO_ROW);
                  removerow[r1] = TRUE;
                  nzerorowsremoved++;
               }
            }
         if( removelargeslackrows && !removerow[r1] )
            if( SCIPisGT(scip, mod2data->slacks[mod2data->rowsind[firstrowsind + r1]], maxslack) )
            {
               markRowAsRemoved(mod2data, firstrowsind + r1, SLACK_GREATER_THAN_MAXSLACK);
               removerow[r1] = TRUE;
               nlargeslackrowsremoved++;
            }

         if( removeidenticalrows && !removerow[r1] )
            for( r2 = r1 + 1 ; r2 < lastrowsind - firstrowsind ; ++r2)
               if( !rowisprocessed[r2] )
                  if( mod2data->rhs[mod2data->rowsind[firstrowsind + r1]]
                     == mod2data->rhs[mod2data->rowsind[firstrowsind + r2]] )
                     if( BITARRAYSAREEQUAL(mod2data->rows[mod2data->rowsind[firstrowsind + r1]],
                           mod2data->rows[mod2data->rowsind[firstrowsind + r2]], mod2data->rowsbitarraysize) )
                     {
                        if( SCIPisLT(scip, mod2data->slacks[mod2data->rowsind[firstrowsind + r1]],            
                              mod2data->slacks[mod2data->rowsind[firstrowsind + r2]]) )
                        {
                           markRowAsRemoved(mod2data, firstrowsind + r2, IDENT_TO_ROW_WITH_SMALLER_SLACK);
                           removerow[r2] = TRUE;
                           nidenticalrowsremoved++;
                           rowisprocessed[r2] = TRUE;
                        }
                        else
                        {
                           markRowAsRemoved(mod2data, firstrowsind + r1, IDENT_TO_ROW_WITH_SMALLER_SLACK);
                           removerow[r1] = TRUE;
                           nidenticalrowsremoved++;
                           break;
                        }  
                     }
      }

   /* update mod2data->rowsind if necessary */
   if( nzerorowsremoved + nlargeslackrowsremoved + nidenticalrowsremoved > 0 )
   {   
      r2 = firstrowsind;
      for( r1 = firstrowsind ; r1 < mod2data->nrowsind && r2 < mod2data->nrowsind; ++r1)
      {
         if( r1 < lastrowsind - firstrowsind )
            while( removerow[r2] && r2 < mod2data->nrowsind )
               r2++;
         if( r1 < r2 && r2 < mod2data->nrowsind )
            mod2data->rowsind[r1] = mod2data->rowsind[r2];
         r2++;
      }      
      mod2data->nrowsind -= (nzerorowsremoved + nlargeslackrowsremoved + nidenticalrowsremoved);
   }
 
  
   /* free temporary memory */
   if( removezerorows && zerorow != NULL )
   {
      SCIPfreeBufferArray(scip, &zerorow);
   }   
   SCIPfreeBufferArray(scip, &removerow);
   SCIPfreeBufferArray(scip, &rowisprocessed);

   return SCIP_OKAY;
}


/** applies some column reductions */
static
SCIP_RETCODE preprocessColumns(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */  
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered (preprocessed) subproblem mod 2 */
   int                   firstcolsind,       /**< first mod2data->rows index to be considered */ 
   int                   lastcolsind,        /**< last mod2data->rows index to be considered */
   SCIP_Bool             removezerocols,     /**< should zero columns be removed? */
   SCIP_Bool             removecolsingletons,/**< should column singletons be removed? */
   SCIP_Bool             checkresultingrows  /**< should rows whose slack becomes larger than maxslack be removed? */
   )
{
   SCIP_Real             maxslack;
   int                   maxnnonzentries;  
   int                   nzerocolsremoved;
   int                   ncolsingletonsremoved;
   int                   nunprocessedcols;
   int                   nconsideredcols;
   int                   nnonzentries;
   SCIP_Bool*            removecol;
   SCIP_Bool*            colisprocessed;  
   int                   rowofcolsingleton;
   int                   rowsbind;
   BITARRAYBITMASKTYPE   rowsbmask;
   int                   r;
   int                   c;
   int                   j;
  
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(lpdata != NULL);
   assert(mod2data != NULL);
   assert(mod2data->relatedsubproblem != NULL);
   assert(firstcolsind >= 0);
   assert(lastcolsind <= mod2data->ncolsind);
   assert(removezerocols || removecolsingletons);
  

   nconsideredcols = lastcolsind - firstcolsind; 
  
   /* check if matrix or colind range is empty */
   if( mod2data->ncolsind == 0 || mod2data->nrowsind == 0 || nconsideredcols <= 0 )
      return SCIP_OKAY;


   /* allocate temporary memory */
   SCIP_CALL(SCIPallocBufferArray(scip, &colisprocessed, nconsideredcols));
   SCIP_CALL(SCIPallocBufferArray(scip, &removecol, nconsideredcols));
  
   /* initialize */
   BMSclearMemoryArray(colisprocessed, nconsideredcols);
   BMSclearMemoryArray(removecol, nconsideredcols);
   maxslack = sepadata->maxslack;
   nunprocessedcols = nconsideredcols;
   nzerocolsremoved = 0;
   ncolsingletonsremoved = 0;
   nnonzentries = 0;
   rowofcolsingleton = -1;
   if( removecolsingletons )
      maxnnonzentries = 1;
   else
      maxnnonzentries = 0;
  
   /* check all columns if they contain exactly one nonzero entry */
   while(nunprocessedcols > 0)
   {
      for( c = 0 ; c < nconsideredcols ; ++c )
         if( colisprocessed[c] == FALSE )
            break;
      assert(firstcolsind + c < mod2data->ncolsind);

      nnonzentries = 0;
      rowsbind = (int) GETBITARRAYINDEX(mod2data->colsind[firstcolsind + c]);
      rowsbmask = GETBITARRAYMASK(mod2data->colsind[firstcolsind + c]); /*lint !e701*/

      for( r = 0 ; r < mod2data->nrowsind ; ++r)
         if( mod2data->rows[mod2data->rowsind[r]][rowsbind] & rowsbmask )
         {
            nnonzentries++;
            if( nnonzentries > maxnnonzentries )
               break;
            rowofcolsingleton = r;
         }

      /* check if a zero column has been found */
      if( removezerocols )
         if( nnonzentries == 0 )
         {
            /* remove zero columns */
            removecol[c] = TRUE;
            nzerocolsremoved++;
            markColAsRemovedAndClearCol(mod2data, firstcolsind + c, ZERO_COLUMN);
         }

      /* check if a column singleton has been found */
      if( removecolsingletons && !removecol[c] )
         if( nnonzentries == 1 )
         {
            r = rowofcolsingleton;
            removecol[c] = TRUE;

            /* update row slack:  slack' = slack + fracsol */
            mod2data->slacks[mod2data->rowsind[r]] += mod2data->fracsol[mod2data->colsind[firstcolsind + c]];

            /* if removing col results in a row with slack > maxslack, 
             * then the row can be removed as well */
            if( checkresultingrows && SCIPisGT(scip, mod2data->slacks[mod2data->rowsind[r]], maxslack) )
            {       
               markRowAsRemoved(mod2data, r, SLACK_GREATER_THAN_MAXSLACK);
               for( j = 0 ; j < nconsideredcols ; ++j) 
                  if( !removecol[j] && colisprocessed[j] )
                     if( BITARRAYBITISSET(mod2data->rows[mod2data->rowsind[r]], 
                           mod2data->colsind[firstcolsind + j]) ) /*lint !e701*/
                     {
                        colisprocessed[j] = FALSE; /* re-consider col */
                        nunprocessedcols++;
                     }

               BMSmoveMemoryArray(&((mod2data->rowsind)[r]), &((mod2data->rowsind)[r + 1]),
                  mod2data->nrowsind - r - 1); /*lint !e866*/

               mod2data->nrowsind--;
            }

            /* remove column singleton */
            ncolsingletonsremoved++;
            markColAsRemovedAndClearCol(mod2data, firstcolsind + c, SINGLETON_COLUMN);
         }
    
      colisprocessed[c] = TRUE;
      nunprocessedcols--;
    
      if( nzerocolsremoved + ncolsingletonsremoved == nconsideredcols || mod2data->nrowsind == 0 )
         break;    
   }

   /* if all rows have been deleted, remove cols as well */
   if( mod2data->nrowsind == 0 )
   {
      for( c = firstcolsind ; c < mod2data->ncolsind; ++c)
         if( !removecol[c] )
         {
            removecol[c] = TRUE;
            markColAsRemovedAndClearCol(mod2data, firstcolsind + c, ZERO_COLUMN);
         }
      nzerocolsremoved = nconsideredcols - ncolsingletonsremoved;
      assert(nzerocolsremoved + ncolsingletonsremoved == nconsideredcols);
   }
    
   /* update mod2data->colsind array if necessary*/
   if( mod2data->nrowsind == 0 )
      mod2data->ncolsind = 0;
   else    
      if( nzerocolsremoved + ncolsingletonsremoved > 0 )
      {    
         j = firstcolsind;
         for( c = firstcolsind ; c < mod2data->ncolsind && j < mod2data->ncolsind; ++c)
         {
            if( c < nconsideredcols )
               while( removecol[j] && j < mod2data->ncolsind )
                  j++;
            if( c < j && j < mod2data->ncolsind )
               mod2data->colsind[c] = mod2data->colsind[j];
            j++;
         }      
         mod2data->ncolsind -= (nzerocolsremoved + ncolsingletonsremoved);
      }
  
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &removecol);
   SCIPfreeBufferArray(scip, &colisprocessed);
  
   return SCIP_OKAY;
}


/** applies modified Gaussian Elimination reduction */
static
SCIP_RETCODE preprocessModGaussElim(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */       
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   ZEROHALF_MOD2DATA*    mod2data            /**< considered (preprocessed) subproblem mod 2 */
   )
{
   int                   nslackzerorows;
   int                   pivotrow;
   int                   pivotcol;
   int                   pivot;
   int                   identsubmatrixsize;
   int                   rowsbind;
   BITARRAYBITMASKTYPE   rowsbmask;
   int                   r;
   int                   temp;
  
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(lpdata != NULL);
   assert(mod2data != NULL);
   assert(mod2data->relatedsubproblem != NULL);

   /* check if matrix or colind range is empty */
   if( mod2data->ncolsind == 0 || mod2data->nrowsind == 0 )
      return SCIP_OKAY;

   /* determine number of slack zero rows */
   nslackzerorows = 0;
   while( nslackzerorows < mod2data->nrowsind
      && SCIPisZero(scip, mod2data->slacks[mod2data->rowsind[nslackzerorows]]) )
      nslackzerorows++;
   /* check if at least one slack zero row exists */
   if( nslackzerorows == 0 )
      return SCIP_OKAY;


   /* sort column indices sets w.r.t. to their primsol values NON-INCREASINGLY */
   if( mod2data->ncolsind > 1 )
   {
      SCIPsortInd( mod2data->colsind , compRealNonIncreasing , (void*) mod2data->fracsol , mod2data->ncolsind );
   }
  
   /* sort row indices sets w.r.t. to their slack values NON-DECREASINGLY */
   if( mod2data->nrowsind > 1 )
   {
      SCIPsortInd( mod2data->rowsind , compRealNonDecreasing , (void*) mod2data->slacks , mod2data->nrowsind );
   }

  
   identsubmatrixsize = 0;
  
   /* create maximal identity submatrix */
   /* determine pivot col */
   for( pivotcol = 0 ; pivotcol < mod2data->ncolsind ; ++pivotcol)
   {
      if( identsubmatrixsize ==  mod2data->nrowsind )
         break;

      /* determine pivot row */
      rowsbind = (int) GETBITARRAYINDEX(mod2data->colsind[pivotcol]);
      rowsbmask = GETBITARRAYMASK(mod2data->colsind[pivotcol]); /*lint !e701*/
      for( pivotrow = identsubmatrixsize ; pivotrow < nslackzerorows ; ++pivotrow)
         if( mod2data->rows[mod2data->rowsind[pivotrow]][rowsbind] & rowsbmask )
            break;
      if( pivotrow == nslackzerorows )
         continue;

      /* Gaussian elimination step */
      for( r = 0 ; r < nslackzerorows ; ++r)
      {
         if( r == pivotrow )
            continue;
         if( mod2data->rows[mod2data->rowsind[r]][rowsbind] & rowsbmask )
         {
            /* add pivot row to r-th row */
            BITARRAYSXOR(mod2data->rows[mod2data->rowsind[pivotrow]],
               mod2data->rows[mod2data->rowsind[r]],mod2data->rowsbitarraysize);
            BITARRAYSXOR(mod2data->rowaggregations[mod2data->rowsind[pivotrow]],
               mod2data->rowaggregations[mod2data->rowsind[r]],mod2data->rowaggregationsbitarraysize);
            mod2data->rhs[mod2data->rowsind[r]] =
               XOR(mod2data->rhs[mod2data->rowsind[pivotrow]],mod2data->rhs[mod2data->rowsind[r]]);
            /* all rows have slack zero: */
            /*            // mod2data->slacks[[mod2data->rowsind[r]] += mod2data->slacks[[mod2data->rowsind[pivotrow]];        */
         }      
      }
      
      /* swap index set positions */
      temp = mod2data->rowsind[pivotrow];
      mod2data->rowsind[pivotrow] = mod2data->rowsind[identsubmatrixsize];
      mod2data->rowsind[identsubmatrixsize] = temp;
      temp = mod2data->colsind[pivotcol];
      mod2data->colsind[pivotcol] = mod2data->colsind[identsubmatrixsize];
      mod2data->colsind[identsubmatrixsize] = temp;

      identsubmatrixsize++;
   }

   if( identsubmatrixsize > 0 )
   {
      /* add rows of identity submatrix properly to rows with positive slack
       * to transform each column of the identity submatrix into a column singleton */
      for( pivot = 0 ; pivot < identsubmatrixsize ; ++pivot)
      {
         rowsbind = (int) GETBITARRAYINDEX(mod2data->colsind[pivot]);
         rowsbmask = GETBITARRAYMASK(mod2data->colsind[pivot]); /*lint !e701*/
         for( r = nslackzerorows ; r < mod2data->nrowsind ; ++r)
            if( mod2data->rows[mod2data->rowsind[r]][rowsbind] & rowsbmask )
            {
               /* add pivot row to r-th row */
               BITARRAYSXOR(mod2data->rows[mod2data->rowsind[pivot]],
                  mod2data->rows[mod2data->rowsind[r]],mod2data->rowsbitarraysize);
               BITARRAYSXOR(mod2data->rowaggregations[mod2data->rowsind[pivot]],
                  mod2data->rowaggregations[mod2data->rowsind[r]],mod2data->rowaggregationsbitarraysize);
               mod2data->rhs[mod2data->rowsind[r]] =
                  XOR(mod2data->rhs[mod2data->rowsind[pivot]],mod2data->rhs[mod2data->rowsind[r]]);
               /* all identity submatrix rows have slack zero */
               /* // mod2data->slacks[[mod2data->rowsind[r]] += mod2data->slacks[[mod2data->rowsind[pivot]]; */
            }              
      }

      /* remove generated column singletons */
      SCIP_CALL(preprocessColumns(scip, sepadata, lpdata, mod2data,
            0, /*identsubmatrixsize*/ mod2data->ncolsind, FALSE, TRUE, TRUE));
        
      /* remove zero rows */
      SCIP_CALL(preprocessRows(scip, sepadata, lpdata, mod2data,
            0, mod2data->nrowsind, TRUE, FALSE, FALSE));
   }
   
   return SCIP_OKAY;
}


/** decomposes the problem into subproblems which can be considered separately */
static
SCIP_RETCODE decomposeProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */       
   ZEROHALF_LPDATA*      lpdata              /**< data of current LP relaxation */  
   )
{

#ifdef WITHDECOMPOSE
   /**@todo this is buggy in different ways.
       * 1. it might happen that we ignore a variable of the current row and of all other rows. 
       * thus at the end, the variable will not occur in any subproblem. BUT, currently we do not update 
       * lpdata->subproblemsindexofcol[lppos] and lpdata->rcolsindexofcol[lppos] accordingly. 
       * consequently, it might happen that lpdata->rcolsindexofcol[lppos] > problem->nrcols, with 
       * with problem being the subproblem still associated to our column. therefore, a corresponding assert 
       * assert(rcolsindex < problem->nrcols) in storeMod2Data() is violated 
       * [e.g., for IP/atamtuerk/mik/unbounded/mik.250-1-50.3.mps.gz].
       * we could recognize whether a variable is never added to a subproblem and update its data structures,
       * but I'm not sure whether this will be correct, i.e., whether it is really ok to ignore some variables here.
       * 2. in particular, it seems like that this method has not been adapted to that we can deal with 
       * continuous variables now. (see the below todo of Manuel)
       * 3. in case we end up with only one subproblem, we use the old problem. but in this case we do not update 
       * the problem data and hence it is not consistent with the lpdata anymore where we might have set some 
       * rows to be irrelevant.  
       *
       * therefore, we will currently do nothing in here.
       */
   BITARRAY              processedrows;
   int                   nprocessedrows;
   int                   processedrowsbitarraysize;
   int                   unprocessedrowidx;

   BITARRAY              processedcols;
   int                   nprocessedcols;
   int                   processedcolsbitarraysize;
  
   int*                  queue;
   int                   queuefirst;
   int                   queuelast;

   int                   i; 
   int                   j;
   int                   k;
  
   SCIP_COL**            colsofrow;
   SCIP_ROW**            rowsofcol;

   SCIP_Real*            colvals;
   SCIP_Real*            rowvals;
   int                   ncolvals;
   int                   nrowvals;
   int                   cidx;
   int                   ridx;
   int                   lppos;
  
   int                   rrowsidx;
   int                   rcolsidx;
  
   SCIP_Bool             fliplhsrhs;
  
   ZEROHALF_SUBLPDATA*   problem;
   int                   problemindex;
   ZEROHALF_SUBLPDATA*   subproblem;
  
   int*                  rrowsinsubprob;
   int*                  rcolsinsubprob;
   SCIP_Bool*            rrowsinsubproboddrhs;
   int                   nrrowsinsubprob;
   int                   nrcolsinsubprob;

   int                   totalnrrows;
   int                   totalnrcols;

   int                   nrrowsinitial;
   int                   nrcolsinitial;
   int                   ndelvarbounds;

   int                   rowindex;    
   int                   colindex;

   SCIP_Real             maxslack;
  
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(lpdata != NULL);
   assert(lpdata->subproblems != NULL);
   assert(lpdata->nsubproblems > 0);

   problemindex = 0;
   
   assert(problemindex >= 0);
   assert(problemindex <= lpdata->nsubproblems);

   problem = lpdata->subproblems[problemindex];
  
   assert(problem != NULL);
   assert(problem->rcols != NULL);
   assert(problem->nrcols > 0);
   assert(problem->rcolslbslack != NULL);
   assert(problem->rcolsubslack != NULL);  
   assert(problem->rrows != NULL);
   assert(problem->nrrows > 0);
   assert(problem->rrowsrhs != NULL);
   assert(problem->rrowsslack != NULL);
  
   if( sepadata->dtimer == NULL )
   {
      ZEROHALFcreateTimer((sepadata->dtimer));    
   }
   ZEROHALFstartTimer(sepadata->dtimer);
    
   processedrowsbitarraysize = (int) GETREQUIREDBITARRAYSIZE(problem->nrrows);
   processedcolsbitarraysize = (int) GETREQUIREDBITARRAYSIZE(problem->nrcols);

   SCIPfreeMemoryArray(scip, &(lpdata->subproblems));

   /* allocate temporary memory */
   SCIP_CALL(SCIPallocMemoryArray(scip, &(lpdata->subproblems), problem->nrrows));  
   SCIP_CALL(SCIPallocBufferArray(scip, &processedrows, processedrowsbitarraysize));
   SCIP_CALL(SCIPallocBufferArray(scip, &processedcols, processedcolsbitarraysize));
   SCIP_CALL(SCIPallocBufferArray(scip, &queue, problem->nrrows));
   SCIP_CALL(SCIPallocMemoryArray(scip, &rrowsinsubprob, problem->nrrows));
   SCIP_CALL(SCIPallocMemoryArray(scip, &rrowsinsubproboddrhs, problem->nrrows));
   SCIP_CALL(SCIPallocMemoryArray(scip, &rcolsinsubprob, problem->nrcols));
  
   /* initialize temporary memory */
   BMSclearMemoryArray(processedrows, processedrowsbitarraysize);
   BMSclearMemoryArray(processedcols, processedcolsbitarraysize);
   lpdata->nsubproblems = 0;
   maxslack = sepadata->maxslack;

   nrrowsinitial = problem->nrrows;
   nrcolsinitial = problem->nrcols;  
   totalnrrows = 0;
   totalnrcols = 0;
   ndelvarbounds = 0;
   nprocessedrows = 0;
   nprocessedcols = 0;
   k = 0;
   unprocessedrowidx = 0;

   while( nprocessedrows < problem->nrrows )
   {
      ++k;
      nrrowsinsubprob = 0;
      nrcolsinsubprob = 0;
    
      for( i = unprocessedrowidx ; i < problem->nrrows ; ++i)
      {
         if( BITARRAYBITISSET(processedrows, i) )
            unprocessedrowidx++;
         else
            break;
      }
      BITARRAYBITSET(processedrows, i);
    
      queue[0] = i;
      queuefirst = 0;
      queuelast = 1;

      while( queuelast > queuefirst )
      {
         assert(queuelast <= problem->nrrows);

         i = queue[queuefirst];
         queuefirst++;

         rrowsinsubprob[nrrowsinsubprob] = i;
         nrrowsinsubprob++;
      
         fliplhsrhs = FALSE;

         colsofrow = SCIProwGetCols(lpdata->rows[problem->rrows[i]]);
         colvals = SCIProwGetVals(lpdata->rows[problem->rrows[i]]);
         ncolvals = SCIProwGetNLPNonz(lpdata->rows[problem->rrows[i]]);

         for( cidx = 0 ; cidx < ncolvals ; ++cidx)
         {
            lppos = SCIPcolGetLPPos(colsofrow[cidx]);
            if( lppos == -1 )
               continue;
            rcolsidx = lpdata->rcolsindexofcol[lppos];

            if( lpdata->subproblemsindexofcol[lppos] != problemindex || rcolsidx < 0 )
            {
               if( ISODD(scip, colvals[cidx]) )          
                  fliplhsrhs = XOR(fliplhsrhs,
                     (lpdata->rcolsindexofcol[lppos] == LP_SOL_EQUALS_ODD_LB
                        || lpdata->rcolsindexofcol[lppos] == LP_SOL_EQUALS_ODD_UB));

               /**@todo analogue for continuous variables? */

               continue;  /* col is not relevant */
            }
            if( nprocessedcols == problem->nrcols )
               continue;
        
            if( BITARRAYBITISSET(processedcols, rcolsidx) )
               continue;

            if( ISEVEN(scip, colvals[cidx]) )
               continue;
        
            rcolsinsubprob[nrcolsinsubprob] = rcolsidx;
            nrcolsinsubprob++;
            BITARRAYBITSET(processedcols, rcolsidx);
        
            rowsofcol = SCIPcolGetRows(colsofrow[cidx]);
            rowvals = SCIPcolGetVals(colsofrow[cidx]);
            nrowvals = SCIPcolGetNNonz(colsofrow[cidx]);
            for( ridx = 0 ; ridx < nrowvals ; ++ridx)
            {
               lppos = SCIProwGetLPPos(rowsofcol[ridx]);
               if( lppos == -1 )
                  continue;
               rrowsidx = lpdata->rrowsindexofleftrow[lppos];
               if( lpdata->subproblemsindexofrow[lppos] == problemindex
                  && rrowsidx >= 0 )
               {
                  if( !BITARRAYBITISSET(processedrows, rrowsidx) )
                     if( ISODD(scip, rowvals[ridx]) )
                     {
                        queue[queuelast] = rrowsidx;
                        queuelast++;
                        BITARRAYBITSET(processedrows, rrowsidx);
                     }
               }
               rrowsidx = lpdata->rrowsindexofrightrow[lppos];
               if(  lpdata->subproblemsindexofrow[lppos] == problemindex
                  && rrowsidx >= 0 )
               {
                  if( !BITARRAYBITISSET(processedrows, rrowsidx) )
                     if( ISODD(scip, rowvals[ridx]) )
                     {
                        queue[queuelast] = rrowsidx;
                        queuelast++;
                        BITARRAYBITSET(processedrows, rrowsidx);
                     }
               }
            }
            nprocessedcols++;
         }
         rrowsinsubproboddrhs[nrrowsinsubprob-1] = XOR(ISODD(scip, problem->rrowsrhs[i]),fliplhsrhs);
         nprocessedrows++;
      }

      /* a subproblem consisting only of rows with even rhs values can be ignored.
       * note: varbounds have to be considered! */
      for( i = 0 ; i < nrrowsinsubprob ; ++i)
         if( rrowsinsubproboddrhs[i] )
            break;
      for( j = 0 ; j < nrcolsinsubprob ; ++j)
         if( (SCIPisLE(scip, problem->rcolsubslack[rcolsinsubprob[j]], maxslack)
               && ISODD(scip, SCIPcolGetUb(lpdata->cols[problem->rcols[rcolsinsubprob[j]]])))
            || (SCIPisLE(scip, problem->rcolslbslack[rcolsinsubprob[j]], maxslack)
               && ISODD(scip, SCIPcolGetLb(lpdata->cols[problem->rcols[rcolsinsubprob[j]]]))) )
            break; /* a relevant odd varbound has been found */    
      if( i == nrrowsinsubprob && j == nrcolsinsubprob )
      {
         /* no odd rhs value exists */
         for( i = 0 ; i < nrrowsinsubprob ; ++i)
         {
            lppos = SCIProwGetLPPos(lpdata->rows[problem->rrows[rrowsinsubprob[i]]]);
            lpdata->subproblemsindexofrow[lppos] = IRRELEVANT;
            if( lpdata->rrowsindexofleftrow[lppos] == rrowsinsubprob[i] )
               lpdata->rrowsindexofleftrow[lppos] = ROW_IN_SUBPROB_WITHOUT_ODD_RHS;
            if( lpdata->rrowsindexofrightrow[lppos] == rrowsinsubprob[i] )
               lpdata->rrowsindexofrightrow[lppos] = ROW_IN_SUBPROB_WITHOUT_ODD_RHS;
         }
         for( j = 0 ; j < nrcolsinsubprob ; ++j)
         {
            lppos = SCIPcolGetLPPos(lpdata->cols[problem->rcols[rcolsinsubprob[j]]]);
            lpdata->subproblemsindexofcol[lppos] = IRRELEVANT;
            lpdata->rcolsindexofcol[lppos] = COLUMN_IN_SUBPROB_WITHOUT_ODD_RHS;                

            /* statistics */
            if( !SCIPisInfinity(scip, problem->rcolslbslack[rcolsinsubprob[j]]) ) 
               ndelvarbounds++;
            if( !SCIPisInfinity(scip, problem->rcolsubslack[rcolsinsubprob[j]]) )
               ndelvarbounds++;
         }
         continue;
      }
    
      /* don't create new "sub"problem if problem can't be decomposed */
      if( lpdata->nsubproblems == 0 && nprocessedrows == problem->nrrows )
         continue; 

      /* create new subproblem */    
      SCIP_CALL(ZerohalfSubLPDataCreate(scip, &subproblem));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(subproblem->rrows), nrrowsinsubprob));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(subproblem->rrowsrhs), nrrowsinsubprob));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(subproblem->rrowsslack), nrrowsinsubprob));
      subproblem->nrrows = nrrowsinsubprob;
      SCIP_CALL(SCIPallocMemoryArray(scip, &(subproblem->rcols), nrcolsinsubprob));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(subproblem->rcolslbslack), nrcolsinsubprob));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(subproblem->rcolsubslack), nrcolsinsubprob));
      subproblem->nrcols = nrcolsinsubprob;

      for( i = 0 ; i < nrrowsinsubprob ; ++i)
      {
         rowindex = problem->rrows[rrowsinsubprob[i]];

         subproblem->rrows[i] = rowindex;
         subproblem->rrowsrhs[i] = problem->rrowsrhs[rrowsinsubprob[i]];    
         subproblem->rrowsslack[i] = problem->rrowsslack[rrowsinsubprob[i]];

         if( lpdata->subproblemsindexofrow[rowindex] != IRRELEVANT )
         {
            assert(lpdata->rrowsindexofleftrow[rowindex] >= 0
               || lpdata->rrowsindexofrightrow[rowindex] >= 0);
            lpdata->subproblemsindexofrow[rowindex] = lpdata->nsubproblems;
            if( lpdata->rrowsindexofleftrow[rowindex] >= 0 )
               lpdata->rrowsindexofleftrow[rowindex] = i;
            if( lpdata->rrowsindexofrightrow[rowindex] >= 0 )
               lpdata->rrowsindexofrightrow[rowindex] = i;
         }
      }

      for( i = 0 ; i < nrcolsinsubprob ; ++i)
      {
         colindex = problem->rcols[rcolsinsubprob[i]];

         subproblem->rcols[i] = colindex;
         subproblem->rcolslbslack[i] = problem->rcolslbslack[rcolsinsubprob[i]];
         subproblem->rcolsubslack[i] = problem->rcolsubslack[rcolsinsubprob[i]];

         if( lpdata->subproblemsindexofcol[colindex] != IRRELEVANT )
         {
            assert(lpdata->rcolsindexofcol[colindex] >= 0);
            lpdata->subproblemsindexofcol[colindex] = lpdata->nsubproblems;        
            lpdata->rcolsindexofcol[colindex] = i;
         }
      }

      lpdata->subproblems[lpdata->nsubproblems] = subproblem;
      lpdata->nsubproblems++;

      totalnrrows += subproblem->nrrows;
      totalnrcols += subproblem->nrcols;

      SCIPdebugMessage("subproblem %d: %d rrows, %d rcols\n", k, subproblem->nrrows, subproblem->nrcols);
   }
   if( lpdata->nsubproblems == 0 )
   {
      /* problem couldn't be decomposed into different subproblems, hence keep the entire problem */
      lpdata->subproblems[0] = problem;
      lpdata->nsubproblems = 1;
      totalnrrows = problem->nrrows;
      totalnrcols = problem->nrcols;
      
   }
   else
   {
      ZerohalfSubLPDataFree(scip, &problem);  
   }
  
   /* free temporary memory */
   SCIPfreeMemoryArray(scip, &rcolsinsubprob);
   SCIPfreeMemoryArray(scip, &rrowsinsubproboddrhs);
   SCIPfreeMemoryArray(scip, &rrowsinsubprob);
   SCIPfreeBufferArray(scip, &queue);
   SCIPfreeBufferArray(scip, &processedcols);
   SCIPfreeBufferArray(scip, &processedrows);


   ZEROHALFstopTimer(sepadata->dtimer);
   ZEROHALFstatisticsMessage("\n");
   ZEROHALFstatisticsMessage("                | --------------------------------- problem \
-------------------------------- | ----- callback ---- | --total-\n");
   ZEROHALFstatisticsMessage("                |   nrrows |   nrcols | ndlrrows | ndlrcols \
| nsubprob | ndelsubp | ndlvbnds | nsepcuts | ncutsfnd |     time\n");
   ZEROHALFstatisticsMessage("%15s | %8d | %8d | %8d | %8d | %8d | %8d | %8d | %8d | %8d | %8.4f\n",
      "DECOMPOSITION", totalnrrows, totalnrcols, nrrowsinitial - totalnrrows, nrcolsinitial - totalnrcols,
      lpdata->nsubproblems, k - lpdata->nsubproblems, 
      ndelvarbounds,
      0, 0, ZEROHALFevalTimer(sepadata->dtimer));
    
#else
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(lpdata != NULL);
#endif

   return SCIP_OKAY;
}


/** removes the largest number of columns such that the sum of the corresponding variables is at most delta */
static
SCIP_RETCODE preprocessColumnsWithSmallFracsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */       
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered (preprocessed) subproblem mod 2 */
   SCIP_Real             delta               /**< delta value */
   )
{
   int                   ncolsremoved;
   int                   c;
   SCIP_Real             maxsumfracsols;
   SCIP_Real             sumfracsols;
  
  
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(mod2data != NULL);
   assert(delta >= 0.0);
   assert(delta <= 1.0);

  
   /* check if matrix contains rows or columns */
   if( mod2data->ncolsind == 0 || mod2data->nrowsind == 0 )
      return SCIP_OKAY;

   /* check if delta is positive */
   if( !SCIPisPositive(scip, delta) )
      return SCIP_OKAY;

    
   ncolsremoved = 0;
   sumfracsols = 0.0;
   maxsumfracsols = sepadata->maxslack * delta;
  
   /* sort column indices sets w.r.t. to their primsol values NON-INCREASINGLY */
   if( mod2data->ncolsind > 1 )
   {
      SCIPsortInd( mod2data->colsind , compRealNonIncreasing , (void*) mod2data->fracsol , mod2data->ncolsind );
   }

   for( c = mod2data->ncolsind - 1 ; c >= 0 ; --c)
   {
      if( SCIPisGT(scip, sumfracsols + mod2data->fracsol[mod2data->colsind[c]], maxsumfracsols) )
         break;

      sumfracsols += mod2data->fracsol[mod2data->colsind[c]];
      markColAsRemovedAndClearCol(mod2data, c, SMALL_FRACSOL_HEUR);
      ncolsremoved++;
   }

   if( ncolsremoved > 0 )
   {  
      mod2data->ncolsind -= ncolsremoved;
      sepadata->maxslack -= sumfracsols;
   }
  
   return SCIP_OKAY;
}


/** removes some rows that cannot be combined because the resulting slack would be larger than maxslack */
static
SCIP_RETCODE preprocessConsiderMinSlack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */  
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered (preprocessed) subproblem mod 2 */ 
   SCIP_Bool             removelargeslackrows, /**< should rows with slack + minslack > maxslack be removed? */
   SCIP_Bool             removelargecolrows  /**< should rows with "large valued" columns that cannot be negated be removed? */
   )
{
   int                   first;
   int                   last;
   int                   temp;
   SCIP_Real             minslackoddrhsrows;
   SCIP_Real             minslackrowwithnonz;
   int                   noddrhsrows;
   int                   nlslrowsremoved;  
   int                   nlcolrowsremoved;
#ifndef NDEBUG
   int                   r;
#endif
   int                   c;
   SCIP_Bool*            removerow;
   int                   i;
   int                   j;
   int                   rowsbind;
   BITARRAYBITMASKTYPE   rowsbmask;
    
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(lpdata != NULL);
   assert(mod2data != NULL);
   assert(removelargeslackrows || removelargecolrows);
  
  
   /* check if( A mod 2, b mod 2) is empty */
   if( mod2data->nrows == 0 || mod2data->nrowsind == 0 )
      return SCIP_OKAY;

  
   /* partition rows into odd-rhs-rows and even-rhs-rows */
   first = 0;
   last = mod2data->nrowsind - 1;
   while(first < last)
   {
      if( !mod2data->rhs[mod2data->rowsind[first]] )
      {
         temp = mod2data->rowsind[first];
         mod2data->rowsind[first] = mod2data->rowsind[last];
         mod2data->rowsind[last] = temp;
         --last;      
      }
      else      
         ++first;                                                                                                             
   } 
   noddrhsrows = first + (mod2data->rhs[mod2data->rowsind[first]] ? 1  : 0);


   /* check if odd rows exists */
   if( noddrhsrows == 0 )
      return SCIP_OKAY;
  
   /* sort each partition by nondecreasing slacks */
   assert(noddrhsrows >= 0);
   SCIPsortInd( mod2data->rowsind , compRealNonDecreasing , (void*) mod2data->slacks , noddrhsrows );
   if( noddrhsrows < mod2data->nrowsind )
   {
      SCIPsortInd( mod2data->rowsind + noddrhsrows , compRealNonDecreasing , (void*) mod2data->slacks , 
         mod2data->nrowsind - noddrhsrows );  
   }
  
   minslackoddrhsrows = mod2data->slacks[mod2data->rowsind[0]];
   nlslrowsremoved = 0;
   nlcolrowsremoved = 0;

   if( SCIPisFeasZero(scip, minslackoddrhsrows) )
      return SCIP_OKAY;
  
   /* check if a zerohalf cut may be generated */
   if( SCIPisGT(scip, minslackoddrhsrows, sepadata->maxslack) )
   {
      if( removelargeslackrows )
      {
         for( i = 0 ; i < mod2data->nrowsind ; ++i )
            markRowAsRemoved(mod2data, i, SLACK_GREATER_THAN_MSL_MINUS_SODD);
         for( i = 0 ; i < mod2data->ncolsind ; ++i )
            markColAsRemovedAndClearCol(mod2data, i, ALL_MATRIX_ROWS_DELETED);
         mod2data->nrowsind = 0;
         mod2data->ncolsind = 0;
      }
   }  
   else
   {    
      SCIP_CALL(SCIPallocBufferArray(scip, &removerow, mod2data->nrowsind));
      BMSclearMemoryArray(removerow,  mod2data->nrowsind);

      /* remove all rows with even rhs and   slack > maxslack - minslackoddrhsrows */
      if( removelargeslackrows )
      {     
         for( i = noddrhsrows ; i < mod2data->nrowsind ; ++i)
            if( SCIPisGT(scip, minslackoddrhsrows + mod2data->slacks[mod2data->rowsind[i]], sepadata->maxslack) )
               break;
         nlslrowsremoved += mod2data->nrowsind - i;      
         while(i < mod2data->nrowsind)
         {
#ifndef NDEBUG
            r = mod2data->rowsind[i];
#endif
            assert(!mod2data->rhs[r]);
            assert(SCIPisGT(scip, minslackoddrhsrows + mod2data->slacks[r], sepadata->maxslack)); 
            markRowAsRemoved(mod2data, i, SLACK_GREATER_THAN_MSL_MINUS_SODD);
            removerow[i] = TRUE;
            i++;
         }    
      }

      /* consider cols */
      if( removelargecolrows )
      {
         if( mod2data->ncolsind > 0 )
         {
            /*   sort column indices sets w.r.t. to their primsol values NON-INCREASINGLY */
            if( mod2data->ncolsind > 1 )
            {
               SCIPsortInd( mod2data->colsind , compRealNonIncreasing , (void*) mod2data->fracsol , mod2data->ncolsind );
            }
        
            j = 0;
            while( j < mod2data->ncolsind && SCIPisGT(scip, mod2data->fracsol[mod2data->colsind[j]] + minslackoddrhsrows, sepadata->maxslack) )
            {
               c = mod2data->colsind[j];
               minslackrowwithnonz = 1.0;
               rowsbind = (int) GETBITARRAYINDEX(c);
               rowsbmask = GETBITARRAYMASK(c); /*lint !e701*/
               for( i = 0 ; i < mod2data->nrowsind ; ++i)
                  if( !removerow[i] )
                     if( mod2data->rows[mod2data->rowsind[i]][rowsbind] & rowsbmask )
                        if( SCIPisLT(scip, mod2data->slacks[mod2data->rowsind[i]], minslackrowwithnonz) )
                           minslackrowwithnonz = mod2data->slacks[mod2data->rowsind[i]];
          
               if( minslackrowwithnonz < 1.0 )
               {
                  for( i = 0 ; i < mod2data->nrowsind ; ++i)
                     if( !removerow[i] )
                        if( mod2data->rows[mod2data->rowsind[i]][rowsbind] & rowsbmask )
                           if( SCIPisGT(scip, minslackrowwithnonz + mod2data->slacks[mod2data->rowsind[i]], sepadata->maxslack) )
                           {
                              markRowAsRemoved(mod2data, i, LARGE_COL_EXISTS);
                              removerow[i] = TRUE;
                              nlcolrowsremoved++;
                           }        
               }        
               j++;
            }      
         }
      }

      /* update mod2data->rowsind if necessary */
      if( nlslrowsremoved + nlcolrowsremoved > 0 )
      {   
         j = 0;
         for( i = 0 ; i < mod2data->nrowsind && j < mod2data->nrowsind; ++i)
         {
            if( i < mod2data->nrowsind )
               while( removerow[j] && j < mod2data->nrowsind )
                  j++;
            if( i < j && j < mod2data->nrowsind )
               mod2data->rowsind[i] = mod2data->rowsind[j];
            j++;
         }      
         mod2data->nrowsind -= (nlslrowsremoved + nlcolrowsremoved);
      }

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &removerow);
   }

   return SCIP_OKAY;
}


/** aggregates identical columns into one column whose (artificial) LP solution is the sum of the aggregated columns */
static
SCIP_RETCODE preprocessIdenticalColums(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_MOD2DATA*    mod2data            /**< considered (preprocessed) subproblem mod 2 */
   )
{
   int                   c1;
   int                   c2;
   int                   r;
   int                   rowsbind1;
   BITARRAYBITMASKTYPE   rowsbmask1;
   int                   rowsbind2;
   BITARRAYBITMASKTYPE   rowsbmask2;
   int                   ncolsremoved;
   SCIP_Bool*            removecol;

   assert(scip != NULL);
   assert(mod2data != NULL);
  
   /* check if( A mod 2, b mod 2) is empty */
   if( mod2data->nrows == 0 || mod2data->nrowsind == 0 || mod2data->ncolsind == 0 )
      return SCIP_OKAY;


   /* allocate and initialize temporary memory */
   SCIP_CALL(SCIPallocBufferArray(scip, &removecol, mod2data->ncolsind));
   BMSclearMemoryArray(removecol, mod2data->ncolsind);
   ncolsremoved = 0;

  
   /* check each pair of columns */
   for( c1 = 0 ; c1 < mod2data->ncolsind - 1 ; ++c1)
   {
      rowsbind1 = (int) GETBITARRAYINDEX(mod2data->colsind[c1]);
      rowsbmask1 = GETBITARRAYMASK(mod2data->colsind[c1]); /*lint !e701*/
      for( c2 = c1 + 1 ; c2 < mod2data->ncolsind ; ++c2)
      {
         rowsbind2 = (int) GETBITARRAYINDEX(mod2data->colsind[c2]);
         rowsbmask2 = GETBITARRAYMASK(mod2data->colsind[c2]); /*lint !e701*/
         for( r = 0 ; r < mod2data->nrowsind ; ++r)
            if( (mod2data->rows[mod2data->rowsind[r]][rowsbind1] & rowsbmask1)
               != (mod2data->rows[mod2data->rowsind[r]][rowsbind2] & rowsbmask2) )
               break;
         if( r == mod2data->nrowsind )
         {
            /* a pair of identical columns have been found */
        
            mod2data->fracsol[mod2data->colsind[c2]] += mod2data->fracsol[mod2data->colsind[c1]];        
            removecol[c1] = TRUE;
            ncolsremoved++;
            markColAsRemovedAndClearCol(mod2data, c1, IDENT_TO_ANOTHER_COLUMN);
         }      
      }
   }

   /* update mod2data->colsind array if necessary*/
   if( ncolsremoved > 0 )
   {    
      c1 = 0;
      for( c2 = 0 ; c2 < mod2data->ncolsind && c1 < mod2data->ncolsind; ++c2)
      {
         if( c2 < mod2data->ncolsind )
            while( removecol[c1] && c1 < mod2data->ncolsind )
               c1++;
         if( c2 < c1 && c1 < mod2data->ncolsind )
            mod2data->colsind[c2] = mod2data->colsind[c1];
         c1++;
      }      
      mod2data->ncolsind -= ncolsremoved;
   }
  
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &removecol);  
  
   return SCIP_OKAY;
}


/** preprocess subproblem */
static
SCIP_RETCODE preprocess(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */   
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered (preprocessed) subproblem mod 2 */
   char                  normtype,           /**< type of norm to use for efficacy norm calculation */
   int                   maxsepacuts,        /**< maximal number of zerohalf cuts separated per separation round */
   int                   maxcuts,            /**< maximal number of zerohalf cuts found per separation round (incl. ineff. cuts) */
   int*                  nsepacuts,          /**< pointer to store current number of separated zerohalf cuts */
   int*                  nzerohalfcuts,      /**< pointer to store current number of found zerohalf cuts */
   ZEROHALF_CUTDATA**    zerohalfcuts,       /**< array to store found zerohalf cuts */
   SCIP_Real**           varsolvals,         /**< dense variable LP solution vector */
   SCIP_RESULT*          result              /**< pointer to SCIP result value of separation */
   )
{
   int                   i;
#ifdef ZEROHALF__PRINT_STATISTICS
   int                   ncolsbeforeppm;
   int                   ncolsinitial;
   int                   nrowsbeforeppm;
   int                   nrowsinitial;
   int                   nsepacutsbeforeppm;
   int                   nsepacutsinitial;
   int                   nzerohalfcutsbeforeppm;
   int                   nzerohalfcutsinitial;
   SCIP_CLOCK* timer;
   SCIP_CLOCK* pptimer;
#endif
   char                  ppname[SCIP_MAXSTRLEN];
  
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(lpdata != NULL);
   assert(mod2data != NULL);
   assert(maxsepacuts >= 0);
   assert(maxcuts >= 0);  
   assert(result != NULL);  
   assert(nsepacuts != NULL);
   assert(nzerohalfcuts != NULL);
   assert(zerohalfcuts != NULL);
   assert(*nsepacuts <= *nzerohalfcuts);
  
   assert(mod2data->relatedsubproblem != NULL);
   assert(mod2data->rows != NULL);
   assert(mod2data->rowaggregations != NULL);
   assert(mod2data->rhs != NULL);
   assert(mod2data->slacks != NULL);
   assert(mod2data->fracsol != NULL);
   assert(mod2data->nrows > 0);
   assert(mod2data->rowsind != NULL);
   assert(mod2data->colsind != NULL);

   if( sepadata->nppmethods == -1 )
   {
      sepadata->nppmethods = (int) strlen(sepadata->ppmethods);
      if( sepadata->nppmethods > 0 && sepadata->ppmethods[0] == '-' )
         sepadata->nppmethods = 0;
   }
  
   if( sepadata->nppmethods == 0 )
      return SCIP_OKAY;

   /* statistics */  
#ifdef ZEROHALF__PRINT_STATISTICS
   if( sepadata->pptimers == NULL )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(sepadata->pptimers), sepadata->nppmethods + 1) );
      for( i = 0 ; i < sepadata->nppmethods + 1 ; ++i)
      {
         ZEROHALFcreateTimer((sepadata->pptimers[i]));
      }
   }
   ZEROHALFstartTimer(sepadata->pptimers[sepadata->nppmethods]);
#endif

   if( mod2data->nrowsind == 0 || mod2data->ncolsind == 0 )
      return SCIP_OKAY;
  
#ifdef ZEROHALF__PRINT_STATISTICS
   ncolsinitial = mod2data->ncolsind;
   nrowsinitial = mod2data->nrowsind;
   nsepacutsinitial = *nsepacuts;
   nzerohalfcutsinitial = *nzerohalfcuts;  
#endif

#ifdef ZEROHALF__PRINT_STATISTICS 
   ZEROHALFstatisticsMessage("\n");
   ZEROHALFstatisticsMessage("                | ------------------------------- subproblem\
 ------------------------------- | ----- callback ---- | --total-\n");
   ZEROHALFstatisticsMessage("                | nrowsind | ncolsind | ndelrows | ndelcols \
| nsepcuts | ncutsfnd |     time | nsepcuts | ncutsfnd |     time\n");
   ZEROHALFstatisticsMessage("%15s | %8d | %8d | %8d | %8d | %8d | %8d | %8.4f | %8d | %8d | %8.4f\n",
      "START PREPROCSS", mod2data->nrowsind, mod2data->ncolsind, 0, 0, 0, 0, 0.0,
      *nsepacuts, *nzerohalfcuts, ZEROHALFevalTimer(sepadata->pptimers[sepadata->nppmethods]));
   ZEROHALFcreateNewTimer(timer);
   ZEROHALFstartTimer(timer);

   ZEROHALFcreateNewTimer(pptimer);  
#endif

   for( i = 0 ; i < sepadata->nppmethods ; ++i)
   {
      /* abort if enough cuts have already been found */
      if( *nsepacuts >= maxsepacuts || *nzerohalfcuts >= maxcuts )
         break;
#ifndef ZEROHALF__PRINT_STATISTICS
      /* abort preprocessing if matrix is empty */
      if( mod2data->nrowsind == 0 && mod2data->ncolsind == 0 )
         break;
#endif

#ifdef ZEROHALF__PRINT_STATISTICS
      /* statistics*/ 
      ZEROHALFstartTimer(pptimer); 
      ZEROHALFstartTimer(sepadata->pptimers[i]);
      ncolsbeforeppm = mod2data->ncolsind;
      nrowsbeforeppm = mod2data->nrowsind;
      nsepacutsbeforeppm = *nsepacuts;
      nzerohalfcutsbeforeppm = *nzerohalfcuts;
#endif
    
      /* apply preprocessing method */    
      switch(sepadata->ppmethods[i])
      {
      case MODGAUSSIANELIMINATION:
         SCIP_CALL(preprocessModGaussElim(scip, sepadata, lpdata, mod2data));
         strncpy(ppname,"gauss",SCIP_MAXSTRLEN);
         break;
      case DELETEZEROROWS:
         SCIP_CALL(preprocessRows(scip, sepadata, lpdata, mod2data,
               0, mod2data->nrowsind, TRUE, FALSE, FALSE));
         strncpy(ppname,"zero rows",SCIP_MAXSTRLEN);
         break;
      case DELETEZEROCOLS:
         SCIP_CALL(preprocessColumns(scip, sepadata, lpdata, mod2data,
               0, mod2data->ncolsind, TRUE, FALSE, FALSE));
         strncpy(ppname,"zero columns",SCIP_MAXSTRLEN);
         break;
      case DELETECOLSINGLETONS:
         SCIP_CALL(preprocessColumns(scip, sepadata, lpdata, mod2data,
               0, mod2data->ncolsind, FALSE, TRUE, TRUE));
         strncpy(ppname,"col singletons",SCIP_MAXSTRLEN);
         break;
      case ADDTRIVIALCUTS:
         SCIP_CALL(preprocessTrivialZerohalfCuts(scip, sepa, sepadata, lpdata, mod2data,
               0, mod2data->nrowsind, normtype, maxsepacuts, maxcuts, nsepacuts,
               nzerohalfcuts, zerohalfcuts, varsolvals, PPZEROONEROW, result));      
         strncpy(ppname,"trivial cuts",SCIP_MAXSTRLEN);
         break;
      case DELETEIDENTROWS:
         SCIP_CALL(preprocessRows(scip, sepadata, lpdata, mod2data,
               0, mod2data->nrowsind, FALSE, FALSE, TRUE));
         strncpy(ppname,"identical rows",SCIP_MAXSTRLEN);
         break;
      case MERGEIDENTCOLS: 
         SCIP_CALL(preprocessIdenticalColums(scip, mod2data));
         strncpy(ppname,"identical cols",SCIP_MAXSTRLEN);
         break;
      case DELETELARGESLACKROWS:
         SCIP_CALL(preprocessRows(scip, sepadata, lpdata, mod2data,
               0, mod2data->nrowsind, FALSE, TRUE, FALSE));
         strncpy(ppname,"sl>maxsl rows",SCIP_MAXSTRLEN);
         break;
      case PPCOLUMNS:
         SCIP_CALL(preprocessColumns(scip, sepadata, lpdata, mod2data,
               0, mod2data->ncolsind, TRUE, TRUE, TRUE));
         strncpy(ppname,"(pp cols)",SCIP_MAXSTRLEN);
         break;
      case PPROWS:
         SCIP_CALL(preprocessRows(scip, sepadata, lpdata, mod2data,
               0, mod2data->nrowsind, TRUE, TRUE, TRUE));
         strncpy(ppname,"(pp rows)",SCIP_MAXSTRLEN);
         break;      
      case DELETESMALLFRACSOLCOLS:
         SCIP_CALL(preprocessColumnsWithSmallFracsol(scip, sepadata, mod2data,
               sepadata->ppdelta));
         strncpy(ppname,"delta heur",SCIP_MAXSTRLEN);
         break;
      case DELETEROWSWRTMINSLACK:
         SCIP_CALL(preprocessConsiderMinSlack(scip, sepadata, lpdata, mod2data, TRUE, TRUE));
         strncpy(ppname,"s_odd",SCIP_MAXSTRLEN);
         break;
      default:
         SCIPerrorMessage("invalid preprocessing method '%c'\n", sepadata->ppmethods[i]);
         return SCIP_INVALIDDATA;   
      }

#ifdef ZEROHALF__PRINT_STATISTICS 
      /* statistics */
      ZEROHALFstopTimer(sepadata->pptimers[i]);
      ZEROHALFstopTimer(pptimer);
      ZEROHALFstatisticsMessage("%15s | %8d | %8d | %8d | %8d | %8d | %8d | %8.4f | %8d | %8d | %8.4f\n",
         ppname, mod2data->nrowsind, mod2data->ncolsind,
         nrowsbeforeppm - mod2data->nrowsind, ncolsbeforeppm - mod2data->ncolsind,
         *nsepacuts - nsepacutsbeforeppm, *nzerohalfcuts - nzerohalfcutsbeforeppm,
         ZEROHALFevalTimer(pptimer), *nsepacuts, *nzerohalfcuts,      
         ZEROHALFevalTimer(sepadata->pptimers[i]));
      ZEROHALFresetTimer(pptimer);    
#endif
   }

#ifdef ZEROHALF__PRINT_STATISTICS 
   /* statistics */
   ZEROHALFstopTimer(timer);
   ZEROHALFstatisticsMessage("%15s | %8d | %8d | %8d | %8d | %8d | %8d | %8.4f | %8d | %8d | %8.4f\n",
      "PREPROCESSED", mod2data->nrowsind, mod2data->ncolsind,
      nrowsinitial - mod2data->nrowsind, ncolsinitial - mod2data->ncolsind,
      *nsepacuts - nsepacutsinitial, *nzerohalfcuts - nzerohalfcutsinitial,
      ZEROHALFevalTimer(timer), *nsepacuts, *nzerohalfcuts,   
      ZEROHALFevalTimer(sepadata->pptimers[sepadata->nppmethods]));
   ZEROHALFstatisticsMessage("\n");
   ZEROHALFfreeTimer(timer);
   ZEROHALFfreeTimer(pptimer);
   ZEROHALFstopTimer(sepadata->pptimers[sepadata->nppmethods]);

   ZEROHALFstatisticsMessage("                | ------------------------------- subproblem ------------------------------- | ------------------------------\n");
   ZEROHALFstatisticsMessage("                |                                           | max2/row | max2/col |  A^T ept |                               \n");
   ZEROHALFstatisticsMessage("%15s |                                           | %8s | %8s | %8s |\n",
      "SUBPROBSTRUCT", 
      hasMatrixMax2EntriesPerRow(mod2data) ? "yes" : "no", hasMatrixMax2EntriesPerColumn(mod2data) ? "yes" : "no", "n/a");
   ZEROHALFstatisticsMessage("\n");
#endif
  
   return SCIP_OKAY;
}


/* --------------------------------------------------------------------------------------------------------------------
 * local methods: separating methods
 * -------------------------------------------------------------------------------------------------------------------- */


/** returns the objective weights for the weighted feasibility AuxIP */
static
SCIP_Real calcObjWeight(
   BITARRAY              rowaggregation,     /**< row aggregation bitarray */
   int                   nrrows              /**< number of relevant rows */
   )
{
   int                   i;
   int                   naggregatedrrows;
  
   assert(rowaggregation != NULL);
   assert(nrrows > 0);
   
   naggregatedrrows = 0;
   for( i = 0 ; i < nrrows ; ++i)
      if( BITARRAYBITISSET(rowaggregation, i) ) /*lint !e701*/
         naggregatedrrows++;
   
   return (SCIP_Real) naggregatedrrows;
}


/** creates a "subscip" representing the following auxiliary IP (AuxIP):
 *  min   z :=           s^T v + x^T y
 *  s.t.       (b (mod 2))^T v         - 2q     = 1 
 *  (A (mod 2))^T v -     y -    -2r = 0 
 *
 *  v \\in {0,1}^nrowsind
 *  y \\in {0,1}^ncolsind
 *  r \\in Z^ncolsind_+
 *  q \\in Z_+  
 */
#define BRANCHPRIORITY__AVOID_BRANCHING    0
#define BRANCHPRIORITY__PREFER_BRANCHING   0
static
SCIP_RETCODE createSubscip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */      
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered (preprocessed) subproblem mod 2 */
   ZEROHALF_AUXIPDATA*   auxipdata,          /**< pointer to data structure to store the auxiliary IP data */
   SCIP_Bool             setnodelimit        /**< should a node limit be set? */
   )
{
   SCIP_VAR**            consvars;
   SCIP_Real*            consvals;
   SCIP_Real             maxslack;
   char                  consname[SCIP_MAXSTRLEN];
   char                  varname[SCIP_MAXSTRLEN];
   int                   i;
   int                   j;
   int                   maxnconsvars;
   int                   nconsvars;
   SCIP_Bool             isfeasip;
   SCIP_Bool             isweighted;
   SCIP_Bool             ispenalized;
   SCIP_Bool             settingsfileexists;

   int                   nrrows;
  
   int                   rowsbind;
   BITARRAYBITMASKTYPE   rowsbmask;

   SCIP_Bool success;

   SCIP_Real             feastol;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(lpdata != NULL);
   assert(mod2data != NULL);
   assert(auxipdata != NULL);
  
   assert(mod2data->relatedsubproblem != NULL);
   assert(mod2data->rows != NULL);
   assert(mod2data->rowaggregations != NULL);
   assert(mod2data->rhs != NULL);
   assert(mod2data->slacks != NULL);
   assert(mod2data->fracsol != NULL);
   assert(mod2data->rowsind != NULL);
   assert(mod2data->colsind != NULL);

   assert(auxipdata->subscip == NULL);
   assert(auxipdata->v == NULL);
   assert(auxipdata->y == NULL);
   assert(auxipdata->r == NULL);
   assert(auxipdata->q == NULL);
   assert(auxipdata->feasipcons == NULL);
   assert(auxipdata->oddrhscons == NULL);
   assert(auxipdata->columnsumcons == NULL);

  
   auxipdata->m = mod2data->nrowsind;
   auxipdata->n = mod2data->ncolsind;
  
   /* alloc temporary memory for subscipdata elements*/
   SCIP_CALL(SCIPallocMemoryArray(scip, &(auxipdata->v), auxipdata->m));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(auxipdata->y), auxipdata->n));  
   SCIP_CALL(SCIPallocMemoryArray(scip, &(auxipdata->r), auxipdata->n));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(auxipdata->columnsumcons), auxipdata->n));

   /* initialize allocated data structures */
   BMSclearMemoryArray(auxipdata->v, auxipdata->m);   /* NULL = 0x0 */
   BMSclearMemoryArray(auxipdata->y, auxipdata->n);   /* NULL = 0x0 */
   BMSclearMemoryArray(auxipdata->r, auxipdata->n);   /* NULL = 0x0 */
   BMSclearMemoryArray(auxipdata->columnsumcons, auxipdata->n);   /* NULL = 0x0 */
  
   maxslack = sepadata->maxslack;
   nrrows = mod2data->relatedsubproblem->nrrows;

   /* determine subscip limits */    
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &auxipdata->timelimit) );
   if( !SCIPisInfinity(scip, auxipdata->timelimit) )
      auxipdata->timelimit -= SCIPgetSolvingTime(scip);

   /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &auxipdata->memorylimit) );
   if( !SCIPisInfinity(scip, auxipdata->memorylimit) )
   {
      auxipdata->memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
      auxipdata->memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
   }

   if( setnodelimit == TRUE )
      auxipdata->nodelimit = 3000;
   else
      auxipdata->nodelimit = -1;

   feastol = SCIPfeastol(scip);
   auxipdata->objectivelimit = MIN(1.0, maxslack + feastol);
  
   /* abort if not enough memory available */
   if( auxipdata->memorylimit <= 2.0*SCIPgetMemExternEstim(scip)/1048576.0 )
      return SCIP_OKAY;

   /* abort if not enough time available */
   if( auxipdata->timelimit <= 0.0 )
      return SCIP_OKAY;

   /* alloc further temporary memory */
   maxnconsvars = auxipdata->m + auxipdata->n + 2;
   SCIP_CALL(SCIPallocBufferArray(scip, &consvals, maxnconsvars));
   SCIP_CALL(SCIPallocBufferArray(scip, &consvars, maxnconsvars));

   /* create and initialize framework */

   SCIP_CALL( SCIPcreate(&(auxipdata->subscip)) ); 
   success = FALSE;
#ifndef NDEBUG
   SCIP_CALL( SCIPcopyPlugins(scip, auxipdata->subscip, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, &success) );
#else
   SCIP_CALL( SCIPcopyPlugins(scip, auxipdata->subscip, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, &success) );
#endif
   SCIPdebugMessage("Copying the plugins was %s successful.\n", success ? "" : "not");

   SCIP_CALL( SCIPcreateProb(auxipdata->subscip, "sepa_zerohalf auxiliary IP (AuxIP)",
         NULL, NULL , NULL , NULL , NULL , NULL , NULL) );

   settingsfileexists = TRUE;
   if( strlen(sepadata->subscipsettings) == 0 )
      settingsfileexists = FALSE;
   if( strlen(sepadata->subscipsettings) == 1 && sepadata->subscipsettings[0] == '-' )
      settingsfileexists = FALSE;

   if( settingsfileexists )
   {
      /* read subscip settings file */
      SCIP_CALL(SCIPreadParams(auxipdata->subscip, sepadata->subscipsettings));
   }
   else
   {
      /* do not abort subscip on CTRL-C */
      SCIP_CALL(SCIPsetBoolParam(auxipdata->subscip, "misc/catchctrlc", FALSE));

      /* disable output to console */
#ifdef SCIP_DEBUG
      SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "display/verblevel", 4));
      SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "display/freq", 1));
#else
      SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "display/verblevel", 0));
      SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "display/freq", 1000));
#endif
      SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "display/nsols/active", 2));

      /* forbid recursive call of heuristics solving subMIPs */
      SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/rins/freq", -1));
      SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/rens/freq", -1));
      SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/localbranching/freq", -1));
      SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/crossover/freq", -1));

      /* disable cut separation in subscip */
      SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "separating/zerohalf/freq", -1));
      /*    SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "separating/maxrounds", 0)); */
      /*     SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "separating/maxroundsroot", 0));  */
      /*     SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "separating/maxcuts", 0));  */
      /*     SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "separating/maxcutsroot", 0)); */ 
    
      /* use pseudo cost branching without strong branching */
      /*     SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "branching/pscost/priority", INT_MAX/4)); */
    
      /* disable expensive presolving */
      /*     SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "presolving/probing/maxrounds", 0)); */
      /*     SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "constraints/linear/maxpresolpairrounds", 0)); */
      /*     SCIP_CALL(SCIPsetRealParam(auxipdata->subscip, "constraints/linear/maxaggrnormscale", 0.0)); */
    
      /* disable conflict analysis */
      /*     SCIP_CALL(SCIPsetBoolParam(auxipdata->subscip, "conflict/useprop", FALSE)); */
      /*     SCIP_CALL(SCIPsetBoolParam(auxipdata->subscip, "conflict/useinflp", FALSE)); */
      /*     SCIP_CALL(SCIPsetBoolParam(auxipdata->subscip, "conflict/useboundlp", FALSE)); */
      /*     SCIP_CALL(SCIPsetBoolParam(auxipdata->subscip, "conflict/usesb", FALSE)); */
      /*     SCIP_CALL(SCIPsetBoolParam(auxipdata->subscip, "conflict/usepseudo", FALSE)); */
    
      SCIP_CALL(SCIPsetBoolParam(auxipdata->subscip, "branching/preferbinary",        TRUE));
      SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/shifting/freq",          3));
      SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/simplerounding/freq",    1));
      SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/rounding/freq",          1));
      SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/oneopt/freq",            1));
    
      /*     SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/pscostdiving/freq",      1)); */
      /*     SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/feaspump/freq",          3)); */
    
      /*     SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/coefdiving/freq",       -1)); */
      /*     SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/fracdiving/freq",       -1)); */
      /*     SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/guideddiving/freq",     -1)); */
      /*     SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/linesearchdiving/freq", -1)); */
      /*     SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/objpscostdiving/freq",  -1)); */
      /*     SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/rootsoldiving/freq",    -1)); */
      /*     SCIP_CALL(SCIPsetIntParam(auxipdata->subscip, "heuristics/veclendiving/freq",     -1)); */
   }
  
   /* get type of auxiliary IP objective function */
   isfeasip = (sepadata->subscipobjective == 'v' ? FALSE : TRUE);
   isweighted = (sepadata->subscipobjective == 'w' ? TRUE : FALSE);
   ispenalized = (sepadata->subscipobjective == 'p' ? TRUE : FALSE);
  
   /* set limits of subscip */
   SCIP_CALL( SCIPsetLongintParam(auxipdata->subscip, "limits/nodes", (SCIP_Longint) auxipdata->nodelimit) );
   SCIP_CALL( SCIPsetRealParam(auxipdata->subscip, "limits/time", auxipdata->timelimit) );
   SCIP_CALL( SCIPsetRealParam(auxipdata->subscip, "limits/memory", auxipdata->memorylimit) );

   if( !isfeasip )
   {
      SCIP_CALL( SCIPsetObjlimit(auxipdata->subscip, auxipdata->objectivelimit) );
   }
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "limits/solutions", sepadata->subscipsollimit) );


   /* create variables and set objective */ 
   /* q */
   SCIP_CALL( SCIPcreateVar(auxipdata->subscip, &(auxipdata->q), "q", 0.0, SCIPinfinity(auxipdata->subscip),
         0.0, SCIP_VARTYPE_INTEGER, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(auxipdata->subscip, auxipdata->q) );
   SCIP_CALL( SCIPchgVarBranchPriority(auxipdata->subscip, auxipdata->q, BRANCHPRIORITY__AVOID_BRANCHING) );
   /* r */
   for( j = 0 ; j < auxipdata->n ; ++j)
   {
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "r_colsind%d(rcols%d)", j, mod2data->colsind[j]);
      SCIP_CALL( SCIPcreateVar(auxipdata->subscip, &(auxipdata->r[j]), varname, 0.0, SCIPinfinity(auxipdata->subscip),
            0.0, SCIP_VARTYPE_INTEGER, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(auxipdata->subscip, auxipdata->r[j]) );
      SCIP_CALL( SCIPchgVarBranchPriority(auxipdata->subscip, auxipdata->q, BRANCHPRIORITY__AVOID_BRANCHING) );
   }
   /* v */
   for( i = 0 ; i < auxipdata->m ; ++i)
   {
      SCIP_Real objcoef;
      assert(mod2data->rows[mod2data->rowsind[i]] != NULL);

      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "v_rowsind%d(rrows%d)", i, mod2data->rowsind[i]);    
      if( isfeasip )
      {
         if( isweighted ) 
         {
            objcoef = calcObjWeight(mod2data->rowaggregations[mod2data->rowsind[i]], nrrows);
         }      
         else if( ispenalized )
         {
            objcoef = mod2data->slacks[mod2data->rowsind[i]]
               + (sepadata->subscipobjpen
                  * calcObjWeight(mod2data->rowaggregations[mod2data->rowsind[i]], nrrows)); 
         }
         else
            objcoef = 1.0;
      }
      else
         objcoef = mod2data->slacks[mod2data->rowsind[i]];
      assert(!SCIPisFeasNegative(scip, objcoef));
      SCIP_CALL( SCIPcreateVar(auxipdata->subscip, &(auxipdata->v[i]), varname, 0.0, 1.0, objcoef,
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(auxipdata->subscip, auxipdata->v[i]) );
      SCIP_CALL(SCIPchgVarBranchPriority(auxipdata->subscip, auxipdata->q, BRANCHPRIORITY__PREFER_BRANCHING));
   }
   /* y */
   for( j = 0 ; j < auxipdata->n ; ++j)
   {
      SCIP_Real objcoef;
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "y_%d(%d)", j, mod2data->colsind[j]);
      if( isfeasip )
         objcoef = 0.0;
      else
         objcoef = mod2data->fracsol[mod2data->colsind[j]];
      SCIP_CALL( SCIPcreateVar(auxipdata->subscip, &(auxipdata->y[j]), varname, 0.0, 1.0, objcoef,
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(auxipdata->subscip, auxipdata->y[j]) );
   }

   /* create constraints */
   /* "feasibility constraint" */ 
   if( isfeasip )
   {
      nconsvars = 0;
      for( i = 0 ; i < auxipdata->m ; ++i)
      {
         consvals[nconsvars] = mod2data->slacks[mod2data->rowsind[i]];
         consvars[nconsvars] = auxipdata->v[i];
         nconsvars++;
      }
      for( j = 0 ; j < auxipdata->n ; ++j)
      {
         consvals[nconsvars] = mod2data->fracsol[mod2data->colsind[j]];
         consvars[nconsvars] = auxipdata->y[j];
         nconsvars++;
      }
      SCIP_CALL( SCIPcreateConsLinear(auxipdata->subscip, &(auxipdata->oddrhscons), "feas",
            nconsvars, consvars, consvals, 0.0, auxipdata->objectivelimit,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(auxipdata->subscip, auxipdata->oddrhscons) );
   }
   /* "odd rhs" */
   nconsvars = 0;
   for( i = 0 ; i < auxipdata->m ; ++i)
      if( mod2data->rhs[mod2data->rowsind[i]] == TRUE )
      {
         consvals[nconsvars] = 1.0;
         consvars[nconsvars] = auxipdata->v[i];
         nconsvars++;
      }
   consvals[nconsvars] = -2.0;
   consvars[nconsvars] = auxipdata->q;
   nconsvars++;
   SCIP_CALL( SCIPcreateConsLinear(auxipdata->subscip, &(auxipdata->oddrhscons), "odd_rhs",
         nconsvars, consvars, consvals, 1.0, 1.0,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(auxipdata->subscip, auxipdata->oddrhscons) );
   /* "column sum" */
   for( j = 0 ; j < auxipdata->n ; ++j)
   {
      nconsvars = 0;
    
      rowsbind = (int) GETBITARRAYINDEX(mod2data->colsind[j]);
      rowsbmask = GETBITARRAYMASK(mod2data->colsind[j]); /*lint !e701*/
      for( i = 0 ; i < auxipdata->m ; ++i) {
         if( mod2data->rows[mod2data->rowsind[i]][rowsbind] & rowsbmask )
         {
            consvals[nconsvars] = 1.0;
            consvars[nconsvars] = auxipdata->v[i];
            nconsvars++;
         }
      }
      consvals[nconsvars] = -1.0;
      consvars[nconsvars] = auxipdata->y[j];
      nconsvars++;
      consvals[nconsvars] = -2.0;
      consvars[nconsvars] = auxipdata->r[j];
      nconsvars++;

      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "col_%d(%d)_sum", j, mod2data->colsind[j]);
      SCIP_CALL( SCIPcreateConsLinear(auxipdata->subscip, &(auxipdata->columnsumcons[j]) , consname, nconsvars, consvars, consvals, 0.0, 0.0,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(auxipdata->subscip, auxipdata->columnsumcons[j]) );
   }
  
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &consvars);
   SCIPfreeBufferArray(scip, &consvals);

   SCIPdebug( SCIP_CALL( SCIPprintOrigProblem(auxipdata->subscip, NULL, NULL, TRUE) ) );
   
   return SCIP_OKAY;   
}


/** solves the auxiliary IP given as subscip */
static
SCIP_RETCODE solveSubscip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */      
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered (preprocessed) subproblem mod 2 */
   ZEROHALF_AUXIPDATA*   auxipdata,          /**< auxiliary IP data */
   SCIP_SOL***           sols,               /**< pointer to store array of solutions */
   int*                  nsols               /**< pointer to store number of solutions */
   )
{
   SCIP_RETCODE          retcode;
   SCIP_STAGE            subscipstage;
   SCIP_Real             maxslack;
   int                   i;
   int                   j;

   assert(sepadata != NULL);
   assert(mod2data != NULL);
   assert(auxipdata != NULL);
   assert(auxipdata->subscip != NULL);
   assert(sols != NULL);
   assert(*sols == NULL);
   assert(*nsols == 0);

   /* solve AuxIP */
   retcode = SCIPsolve(auxipdata->subscip);

   /* errors in solving the subproblem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop. */
   if ( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage(scip, "Error while solving subproblem in zerohalf separator; sub-SCIP terminated with code <%d>\n", retcode);
      *nsols = 0;
      return SCIP_OKAY;
   }

   /* print statistic */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(auxipdata->subscip, NULL) ) );

   maxslack = sepadata->maxslack;

   /* check if solving was successful and get solutions */
   subscipstage = SCIPgetStage(auxipdata->subscip);
   if( subscipstage == SCIP_STAGE_SOLVING || subscipstage == SCIP_STAGE_SOLVED )
      *nsols = SCIPgetNSols(auxipdata->subscip);
   else
      *nsols = 0;
   if( *nsols > 0 )
   {
      *sols = SCIPgetSols(auxipdata->subscip);
      /* check if only the best solution should be used */
      if( !sepadata->subscipuseallsols )
         *nsols = 1;
   }
  
   /* check if proper a proper solution was found */
   if( sepadata->subscipobjective == 'v' )
   {
      for( i = 0 ; i < *nsols ; ++i)
         if( SCIPisGT(scip, SCIPgetSolOrigObj(auxipdata->subscip , (*sols)[i]), maxslack) )
            break;
      *nsols = i;
   }
   else
   {
      SCIP_Real z;
      SCIP_Real* viols;    
      SCIP_SOL** propersols;

      SCIP_Bool swapped;    
      SCIP_Real tempviol;
      SCIP_SOL* tempsol;

      int npropersols;

      SCIP_CALL(SCIPallocBufferArray(scip, &viols, *nsols));
      SCIP_CALL(SCIPallocMemoryArray(scip, &propersols, *nsols));
      npropersols = 0;
      for( i = 0 ; i < *nsols ; ++i)
      {
         z = 0.0;
         for( j = 0 ; j < auxipdata->m ; ++j)
            z += mod2data->slacks[mod2data->rowsind[j]]
               * SCIPgetSolVal(auxipdata->subscip, (*sols)[i], auxipdata->v[j]);
         for( j = 0 ; j < auxipdata->n ; ++j)
            z += mod2data->fracsol[mod2data->colsind[j]]
               * SCIPgetSolVal(auxipdata->subscip, (*sols)[i], auxipdata->y[j]);

         if( SCIPisLE(scip, z, maxslack) )
         {
            /* proper sol has been found */
            propersols[npropersols] = (*sols)[i];
            viols[npropersols] = z;        
            npropersols++;
         }
      }

    
      swapped = TRUE;
      for( i = 1 ; i < npropersols && swapped; ++i)
      {
         swapped = FALSE;
         for( j = 0 ; j < npropersols - i ; ++j)
         {
            if( viols[j] > viols[j+1] )
            {
               tempviol = viols[j+1];
               viols[j+1] = viols[j];
               viols[j] = tempviol;
               tempsol = propersols[j+1];
               propersols[j+1] = propersols[j];
               propersols[j] = tempsol;            
               swapped = TRUE;
            }
         }
      }
    
      *sols = propersols;
      *nsols = npropersols;
      SCIPfreeBufferArray(scip, &viols);
      
   }
 
   return SCIP_OKAY;
}


  

/** determines the weightvector for a single row */
static
SCIP_RETCODE getZerohalfWeightvectorForSingleRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   int                   rowsindex,          /**< lpdata->rows index */
   int                   rrowsindex,         /**< "subproblem"->rrows index */
   SCIP_Real**           weights             /**< pointer to store the weight vector */
   )
{   /*lint --e{438}*/
   
   assert(scip != NULL);
   assert(lpdata != NULL);
   assert(lpdata->nrows > 0);
   assert(0 <= rowsindex);
   assert(rowsindex < lpdata->nrows);  
   assert(lpdata->rrowsindexofleftrow[rowsindex] == rrowsindex
      || lpdata->rrowsindexofrightrow[rowsindex] == rrowsindex);
   assert(weights != NULL);
   assert(*weights == NULL);


   /* allocate temporary memory */ 
   SCIP_CALL(SCIPallocMemoryArray(scip, weights, lpdata->nrows));

   /* initialize */
   BMSclearMemoryArray(*weights, lpdata->nrows);

   /* determine row weights */
   if( lpdata->rrowsindexofleftrow[rowsindex] == rrowsindex )
      (*weights)[rowsindex] = lpdata->intscalarsleftrow[rowsindex] * (-0.5);
   else
      (*weights)[rowsindex] = lpdata->intscalarsrightrow[rowsindex] * 0.5;

   if( SCIProwGetNLPNonz(lpdata->rows[rowsindex]) >= sepadata->maxnnonz )
   {
      SCIPfreeMemoryArray(scip, weights);
      weights = NULL; 
   }

   return SCIP_OKAY;
}


/** gets the subset of rows that should be combined to a violated zerohalf cut */
static
SCIP_RETCODE getBitarrayOfSelectedRows(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered (preprocessed) subproblem mod 2 */
   ZEROHALF_AUXIPDATA*   auxipdata,          /**< auxiliary IP data */
   SCIP_SOL*             solution,           /**< considered solution */
   BITARRAY*             rrowsincut,         /**< pointer to store the subset of rows as bitarray (length: number of relevant rows) */
   int*                  nrrowsincut         /**< number of combined relevant rows */
   )
{
   int                   i;

   assert(scip != NULL);
   assert(mod2data != NULL);
   assert(mod2data->nrowsind > 0);
   assert(auxipdata != NULL);
   assert(auxipdata->subscip != NULL);
   assert(auxipdata->v != NULL);
   assert(solution != NULL);
   assert(rrowsincut != NULL);
   assert(*rrowsincut == NULL);
   assert(nrrowsincut != NULL);

  
   /* allocate and initialize temporary memory for calculating the symmetric difference */
   SCIP_CALL(SCIPallocMemoryArray(scip, rrowsincut, mod2data->rowaggregationsbitarraysize));
   BITARRAYCLEAR(*rrowsincut, mod2data->rowaggregationsbitarraysize);

   *nrrowsincut = 0;
  
   /* calculate symmetric difference of rrowsincut and specific rowaggregations */
   for( i = 0 ; i < mod2data->nrowsind ; ++i)
      if( auxipdata->v[i] != NULL )     
         if( !SCIPisZero(scip, SCIPgetSolVal(auxipdata->subscip, solution, auxipdata->v[i])) )
         {
            BITARRAYSXOR(mod2data->rowaggregations[mod2data->rowsind[i]], (*rrowsincut) ,
               mod2data->rowaggregationsbitarraysize);
            (*nrrowsincut)++;
         }

   return SCIP_OKAY;
}


/** separates violated zerohalf cuts by solving an auxiliary IP. (exact method; exponential time) */
static
SCIP_RETCODE separateBySolvingAuxIP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */        
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered (preprocessed) subproblem mod 2 */
   char                  normtype,           /**< SCIP normtype */
   int                   maxsepacuts,        /**< maximal number of zerohalf cuts separated per separation round */
   int                   maxcuts,            /**< maximal number of zerohalf cuts found per separation round (incl. ineff. cuts) */
   SCIP_Bool             setnodelimit,       /**< should a node limit be set for solving the auxiliary IP? */
   int*                  nsepacuts,          /**< pointer to store current number of separated zerohalf cuts */
   int*                  nzerohalfcuts,      /**< pointer to store current number of found zerohalf cuts */
   ZEROHALF_CUTDATA**    zerohalfcuts,       /**< array to store found zerohalf cuts */
   SCIP_Real**           varsolvals,         /**< dense variable LP solution vector */
   SCIP_RESULT*          result              /**< pointer to SCIP result value of separation */
   )
{ 
   ZEROHALF_AUXIPDATA*   auxipdata;
   SCIP_SOL**            sols;
   int                   nsols;
   int                   s;
   BITARRAY              rrowsincut;
   int                   nrrowsincut;  
   SCIP_Real*            weights;
   int                   nrowsincut;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(lpdata != NULL);
   assert(mod2data != NULL);
   assert(maxsepacuts >= 0);
   assert(maxcuts >= 0);  
   assert(nsepacuts != NULL);
   assert(nzerohalfcuts != NULL);
   assert(zerohalfcuts != NULL);
   assert(*nsepacuts <= *nzerohalfcuts);
   assert(varsolvals != NULL);
   assert(result != NULL);

   assert(mod2data->relatedsubproblem != NULL);
   assert(mod2data->rows != NULL);
   assert(mod2data->rowaggregations != NULL);
   assert(mod2data->rhs != NULL);
   assert(mod2data->slacks != NULL);
   assert(mod2data->fracsol != NULL);
   assert(mod2data->rowsind != NULL);
   assert(mod2data->colsind != NULL);


   /* check if( A mod 2, b mod 2) is empty */
   if( mod2data->nrows == 0 || mod2data->nrowsind == 0 )
      return SCIP_OKAY;

   /* check if enough cuts have been found */
   if( *nsepacuts >= maxsepacuts || *nzerohalfcuts >= maxcuts )
      return SCIP_OKAY;  
  
   /* allocate temporary memory for subscip data structure */
   SCIP_CALL(ZerohalfAuxIPDataCreate(scip, &auxipdata));
  
   /* create subscip */
   SCIP_CALL(createSubscip(scip, sepadata, lpdata, mod2data, auxipdata, setnodelimit));

   /* abort if subscip was not created */
   if( auxipdata->subscip == NULL )
   {
      SCIP_CALL(ZerohalfAuxIPDataFree(scip, &auxipdata));
      return SCIP_OKAY;
   }
  
   /* solve subscip and get solutions yielding a zerohalf cut with violation >= minviolation */
   sols = NULL;
   nsols = 0;
   SCIP_CALL(solveSubscip(scip, sepadata, mod2data, auxipdata, &sols, &nsols));

  
   /* process solutions */
   for( s = 0; s < nsols ; ++s)
   {
      /* check if enough cuts have been found */
      if( *nsepacuts >= maxsepacuts || *nzerohalfcuts >= maxcuts )
         break;    

      /* determine rrows of the related subproblem that have to be combined */
      rrowsincut = NULL;
      SCIP_CALL(getBitarrayOfSelectedRows(scip, mod2data, auxipdata, sols[s],
            &rrowsincut, &nrrowsincut));
      assert(nrrowsincut > 0);
    
      /* calculate rows zerohalf weightvector */
      weights = NULL;
      SCIP_CALL(getZerohalfWeightvectorFromSelectedRowsBitarray(scip, sepadata, lpdata,
            mod2data, rrowsincut, &weights, &nrowsincut));
      if ( weights == NULL )
      {
         SCIPfreeMemoryArray(scip, &rrowsincut);
         continue;
      }
      assert(nrowsincut > 0);


#ifdef SCIP_DEBUG
      SCIP_CALL( debugPrintLPRowsAndCols(scip, lpdata) );
      SCIPdebugMessage("\n");
      debugPrintSubLpData(scip, lpdata, mod2data->relatedsubproblem);
      debugPrintMod2Data(scip, lpdata, mod2data, TRUE);
      SCIPdebugMessage("\n");
      SCIP_CALL( SCIPprintOrigProblem(auxipdata->subscip, NULL, NULL, TRUE) );
      SCIPdebugMessage("\n");
      SCIP_CALL( SCIPprintBestSol(auxipdata->subscip, NULL , FALSE) );
#endif

      /* create zerohalf cut */
      SCIP_CALL(ZerohalfCutDataCreate(scip, &(zerohalfcuts[*nzerohalfcuts]),
            mod2data->relatedsubproblem, mod2data, nrrowsincut, nrowsincut, AUXIP));
      SCIP_CALL(createZerohalfCutFromZerohalfWeightvector(scip, sepa, sepadata,
            lpdata, weights, normtype, *nzerohalfcuts, varsolvals, zerohalfcuts[*nzerohalfcuts]));

    
      /* add cut */
      SCIP_CALL( addZerohalfCutToLP(scip, sepadata, zerohalfcuts[*nzerohalfcuts], nsepacuts, result) );
      (*nzerohalfcuts)++;

      /* free temporary memory */
      SCIPfreeMemoryArray(scip, &weights);
      SCIPfreeMemoryArray(scip, &rrowsincut);
   }
  

   /* free temporary memory */
   if( sepadata->subscipobjective != 'v' )
   {
      SCIPfreeMemoryArray(scip, &sols);
   }
   SCIP_CALL(ZerohalfAuxIPDataFree(scip, &auxipdata));

   return SCIP_OKAY;
}


/** calculates the inner product of mod2data->row and the LP solution */
static
SCIP_RETCODE calcInnerProductOfRowAndFracsol(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered (preprocessed) subproblem mod 2 */
   BITARRAY              row,                /**< considered mod2data->row */
   SCIP_Real             maxinnerproduct,    /**< calculation is aborted if innerproduct >= maxinnerproduct */
   SCIP_Real*            innerproduct        /**< pointer to store the inner product */
   )
{
   int                   c;
   int                   rcolindex;
  
   assert(scip != NULL);
   assert(mod2data != NULL);
   assert(row != NULL);
   assert(maxinnerproduct >= 0);
   assert(innerproduct != NULL);

  
   *innerproduct = 0.0;

   /* check if( A mod 2, b mod 2) is empty */
   if( mod2data->nrows == 0 || mod2data->nrowsind == 0
      || mod2data->ncolsind == 0 )
      return SCIP_OKAY;

   /* calculate the inner product of rows[rowsindex] and fracsol */
   for( c = 0 ; c < mod2data->ncolsind ; ++c)
   {
      rcolindex = mod2data->colsind[c];
      if( BITARRAYBITISSET(row, rcolindex) ) /*lint !e701*/
         *innerproduct += mod2data->fracsol[rcolindex];
      if( SCIPisGT(scip, *innerproduct, maxinnerproduct) )
         break;
   }
  
   return SCIP_OKAY;
}


/** separate violated zerohalf cuts by enumerating possible row combinations. (heuristic; polynomial time) */
static
SCIP_RETCODE separateByEnumerationHeuristics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */        
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered (preprocessed) subproblem mod 2 */
   char                  normtype,           /**< SCIP normtype */
   int                   maxsepacuts,        /**< maximal number of zerohalf cuts separated per separation round */
   int                   maxcuts,            /**< maximal number of zerohalf cuts found per separation round (incl. ineff. cuts) */
   int*                  nsepacuts,          /**< pointer to store current number of separated zerohalf cuts */
   int*                  nzerohalfcuts,      /**< pointer to store current number of found zerohalf cuts */ 
   ZEROHALF_CUTDATA**    zerohalfcuts,       /**< array to store found zerohalf cuts */
   SCIP_Real**           varsolvals,         /**< dense variable LP solution vector */
   SCIP_RESULT*          result,             /**< pointer to SCIP result value of separation */
   int                   maxncombinedrows    /**< maximal number of combined rows; currently only 1 or 2 is supported */
   )
{
   int                   first;
   int                   last;
   int                   temp;
   SCIP_Real             minslackoddrhsrows;
   int                   ncombinedrows;
   int                   noddrhsrows;
   int                   r1;
   int                   r2;
   BITARRAY              combinedrow;
   BITARRAY              rrowsincut;
   SCIP_Real             roundingdownweakening;
   SCIP_Real             slack1;
   SCIP_Real             slack2;
   int                   i;
   int                   j;
   SCIP_Real*            weights;
   int                   nrowsincut;
  
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(lpdata != NULL);
   assert(mod2data != NULL);
   assert(maxsepacuts >= 0);
   assert(maxcuts >= 0);  
   assert(nsepacuts != NULL);
   assert(nzerohalfcuts != NULL);
   assert(zerohalfcuts != NULL);
   assert(*nsepacuts <= *nzerohalfcuts);
   assert(varsolvals != NULL);
   assert(result != NULL);
   assert(maxncombinedrows >= 1);
  
   assert(mod2data->relatedsubproblem != NULL);
   assert(mod2data->rows != NULL);
   assert(mod2data->rowaggregations != NULL);
   assert(mod2data->rhs != NULL);
   assert(mod2data->slacks != NULL);
   assert(mod2data->fracsol != NULL);
   assert(mod2data->rowsind != NULL);
   assert(mod2data->colsind != NULL);


   /* check if( A mod 2, b mod 2) is empty */
   if( mod2data->nrows == 0 || mod2data->nrowsind == 0 )
      return SCIP_OKAY;

   /* check if enough cuts have been found */
   if( *nsepacuts >= maxsepacuts || *nzerohalfcuts >= maxcuts )
      return SCIP_OKAY; 


   /* partition rows into odd-rhs-rows and even-rhs-rows */
   first = 0;
   last = mod2data->nrowsind - 1;
   while(first < last)
   {
      if( !mod2data->rhs[mod2data->rowsind[first]] )
      {
         temp = mod2data->rowsind[first];
         mod2data->rowsind[first] = mod2data->rowsind[last];
         mod2data->rowsind[last] = temp;
         --last;      
      }
      else      
         ++first;
   }
   noddrhsrows = first + (mod2data->rhs[mod2data->rowsind[first]] ? 1  : 0);
  
   /* check if odd rows exists */
   if( noddrhsrows == 0 )
      return SCIP_OKAY; 

   /* allocate and initialize temporary memory for calculating the symmetric difference */
   combinedrow = NULL;
   rrowsincut = NULL;
   if( maxncombinedrows > 1 )
   {
      SCIP_CALL(SCIPallocBufferArray(scip, &combinedrow, mod2data->rowsbitarraysize));
      SCIP_CALL(SCIPallocBufferArray(scip, &rrowsincut, mod2data->rowaggregationsbitarraysize));
      BITARRAYCLEAR(rrowsincut, mod2data->rowaggregationsbitarraysize);
   }  

   /* sort each partition by nondecreasing slacks */
   assert(noddrhsrows >= 0);
   SCIPsortInd( mod2data->rowsind , compRealNonDecreasing , (void*) mod2data->slacks , noddrhsrows );
   
   if( noddrhsrows < mod2data->nrowsind )
   {
      SCIPsortInd( mod2data->rowsind + noddrhsrows , compRealNonDecreasing , (void*) mod2data->slacks , mod2data->nrowsind - noddrhsrows );
   }

   minslackoddrhsrows = mod2data->slacks[mod2data->rowsind[0]];
  
   if( SCIPisLE(scip, minslackoddrhsrows, sepadata->maxslack) )
      for( ncombinedrows = 1 ; ncombinedrows <= maxncombinedrows ; ++ncombinedrows )
      {
         switch( ncombinedrows )
         {
         case 1:
            /* check all rows r1 with rhs(r1) odd */
            for( i = 0 ; i < noddrhsrows ; ++i)
            {
               r1 = mod2data->rowsind[i];
               assert(mod2data->rhs[r1]);
               if( *nzerohalfcuts == maxcuts || *nsepacuts == maxsepacuts )
                  break;        
               slack1 = mod2data->slacks[r1];
               if( SCIPisGT(scip, slack1, sepadata->maxslack) )
                  break; /* because rowsind is sorted */
               SCIP_CALL(calcInnerProductOfRowAndFracsol(scip, mod2data, mod2data->rows[r1],
                     sepadata->maxslack, &roundingdownweakening));
               if( SCIPisLE(scip, roundingdownweakening + slack1, sepadata->maxslack) )
               {
                  /* a violated zerohalf cut has been found */
            
                  /* calculate rows zerohalf weightvector */
                  weights = NULL;
                  SCIP_CALL(getZerohalfWeightvectorFromSelectedRowsBitarray(scip, sepadata, lpdata,
                        mod2data, mod2data->rowaggregations[r1], &weights, &nrowsincut));
                  if( weights == NULL )
                  {
                     continue;
                  }
                  assert(nrowsincut > 0);

                  /* create zerohalf cut */
                  SCIP_CALL(ZerohalfCutDataCreate(scip, &(zerohalfcuts[*nzerohalfcuts]),
                        mod2data->relatedsubproblem, mod2data, 2, nrowsincut, HEURISTICSENUM));
                  SCIP_CALL(createZerohalfCutFromZerohalfWeightvector(scip, sepa, sepadata,
                        lpdata, weights, normtype, *nzerohalfcuts, varsolvals, zerohalfcuts[*nzerohalfcuts]));
            
                  /* add cut */
                  SCIP_CALL( addZerohalfCutToLP(scip, sepadata, zerohalfcuts[*nzerohalfcuts], nsepacuts, result) );
                  (*nzerohalfcuts)++;
            
                  /* free temporary memory */
                  SCIPfreeMemoryArray(scip, &weights);          
               }        
            }      
            break;
         case 2:        
            assert(combinedrow != NULL);
            assert(rrowsincut != NULL);
            if( noddrhsrows == mod2data->nrowsind )
               break;
            if( mod2data->nrowsind < 2 )
               break;
        
            /* check all pairs (r1,r2) with rhs(r1) odd and rhs(r1) even */
            for( i = 0 ; i < noddrhsrows ; ++i)
            {
               r1 = mod2data->rowsind[i];
               assert(mod2data->rhs[r1]);
               if( *nzerohalfcuts == maxcuts || *nsepacuts == maxsepacuts )
                  break;
               slack1 = mod2data->slacks[r1];
               if( SCIPisGT(scip, slack1, sepadata->maxslack) )
                  break; /* because rowsind_odd is sorted */
          
               for( j = noddrhsrows ; j < mod2data->ncolsind ; ++j)
               {
                  r2 = mod2data->rowsind[j];
                  assert(!mod2data->rhs[r2]);
                  if( *nzerohalfcuts == maxcuts || *nsepacuts == maxsepacuts )
                     break;
                  slack2 = mod2data->slacks[r2];          
                  if( SCIPisGT(scip, slack1 + slack2, sepadata->maxslack) )
                     break; /* because rowsind_even is sorted */
                  BMScopyMemoryArray(combinedrow, mod2data->rows[r1], mod2data->rowsbitarraysize);
                  BITARRAYSXOR(mod2data->rows[r2], combinedrow, mod2data->rowsbitarraysize);
                  SCIP_CALL(calcInnerProductOfRowAndFracsol(scip, mod2data, combinedrow,
                        sepadata->maxslack, &roundingdownweakening));
                  if( SCIPisLE(scip, roundingdownweakening + slack1 + slack2, sepadata->maxslack) )
                  {
                     /* a violated zerohalf cut has been found */
              
                     /* determine rrows of the related subproblem that have to be combined */
                     BMScopyMemoryArray(rrowsincut, mod2data->rowaggregations[r1],
                        mod2data->rowaggregationsbitarraysize);
                     BITARRAYSXOR(mod2data->rowaggregations[r2], rrowsincut, 
                        mod2data->rowaggregationsbitarraysize);
              
                     /* calculate rows zerohalf weightvector */
                     weights = NULL;
                     SCIP_CALL(getZerohalfWeightvectorFromSelectedRowsBitarray(scip, sepadata, lpdata,
                           mod2data, rrowsincut, &weights, &nrowsincut));
                     if ( weights == NULL )
                     {
                        continue;
                     }
                     assert(nrowsincut > 0);

                     /* create zerohalf cut */
                     SCIP_CALL(ZerohalfCutDataCreate(scip, &(zerohalfcuts[*nzerohalfcuts]),
                           mod2data->relatedsubproblem, mod2data, 2, nrowsincut, HEURISTICSENUM));
                     SCIP_CALL(createZerohalfCutFromZerohalfWeightvector(scip, sepa, sepadata,
                           lpdata, weights, normtype, *nzerohalfcuts, varsolvals, zerohalfcuts[*nzerohalfcuts]));

                     /* add cut */
                     SCIP_CALL( addZerohalfCutToLP(scip, sepadata, zerohalfcuts[*nzerohalfcuts], nsepacuts, result) );
                     (*nzerohalfcuts)++;

                     /* free temporary memory */
                     SCIPfreeMemoryArray(scip, &weights);
                  }
               }
            }
            break;
         default:
            SCIPerrorMessage("invalid ncombinedrows '%d'\n", ncombinedrows);
            return SCIP_INVALIDDATA;
         }
      }

   /* free temporary memory */
   if( rrowsincut != NULL )
   {
      SCIPfreeBufferArray(scip, &rrowsincut);  
   }   
   if( combinedrow != NULL )
   {    
      SCIPfreeBufferArray(scip, &combinedrow);  
   }
   return SCIP_OKAY;
}


#if 0
/** prints a node of the auxiliary graph */
static
void debugPrintAuxGraphNode(
   ZEROHALF_AUXGRAPH_NODE* node              /**< node to be printed */
   )
{
   int i;

   assert(node != NULL);

   SCIPdebugMessage("\nnode: %p\n", node);
   for( i = 0 ; i < node->nneighbors ; ++i)
   {
      SCIPdebugMessage("  neighbor %4d: %p  weight: %6f  rrow: %4d\n",
         i, node->neighbors[i], node->edgeweights[i], node->relatedrows[i]);
   }
   SCIPdebugMessage("  nneighbors: %d  distance: %6f  previous: %p\n",
      node->nneighbors, node->distance, node->previous);  
}
#endif  


/** adds an edge (and its "copy" w.r.t. the node copies) to the auxiliary graph */
static
SCIP_RETCODE addEdgeToAuxGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_AUXGRAPH*    graph,              /**< auxiliary graph */
   int                   node1index,         /**< start node of edge */
   int                   node2index,         /**< end node of edge */
   SCIP_Bool             isodd,              /**< is the rhs value of the corresponding mod2data->row odd? */
   SCIP_Real             weight,             /**< weight of the edge */
   int                   relatedrow          /**< corresponding mod2data->row */
   )
{
   ZEROHALF_AUXGRAPH_NODE* node1;
   ZEROHALF_AUXGRAPH_NODE* node2;
   ZEROHALF_AUXGRAPH_NODE* node1copy;  
   ZEROHALF_AUXGRAPH_NODE* node2copy;
   int                     n1;
   int                     n2;

   int maxnumberofneighbors;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(node1index >= 0);
   assert(node1index < graph->nnodes);
   assert(node2index >= 0);
   assert(node2index < graph->nnodes);
   assert(!SCIPisNegative(scip, weight));

   maxnumberofneighbors = 2 * graph->nnodes - 2;

   if( isodd )
   {
      node1 = graph->nodes[node1index];
      node2 = graph->nodecopies[node2index];
      node1copy = graph->nodecopies[node1index];
      node2copy = graph->nodes[node2index];
   }
   else
   {
      node1 = graph->nodes[node1index];
      node2 = graph->nodes[node2index];
      node1copy = graph->nodecopies[node1index];
      node2copy = graph->nodecopies[node2index];
   }

   if( node1->nneighbors == 0 )
   {
      assert(maxnumberofneighbors > 0);
      SCIP_CALL(SCIPallocMemoryArray(scip, &(node1->neighbors), maxnumberofneighbors));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(node1->edgeweights), maxnumberofneighbors));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(node1->relatedrows), maxnumberofneighbors));

      SCIP_CALL(SCIPallocMemoryArray(scip, &(node1copy->neighbors), maxnumberofneighbors));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(node1copy->edgeweights), maxnumberofneighbors));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(node1copy->relatedrows), maxnumberofneighbors));
   }

   if( node2->nneighbors == 0 )
   {
      assert(maxnumberofneighbors > 0);
      SCIP_CALL(SCIPallocMemoryArray(scip, &(node2->neighbors), maxnumberofneighbors));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(node2->edgeweights), maxnumberofneighbors));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(node2->relatedrows), maxnumberofneighbors));

      SCIP_CALL(SCIPallocMemoryArray(scip, &(node2copy->neighbors), maxnumberofneighbors));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(node2copy->edgeweights), maxnumberofneighbors));
      SCIP_CALL(SCIPallocMemoryArray(scip, &(node2copy->relatedrows), maxnumberofneighbors));
   }

   n2 = node2->nneighbors;
   for( n1 = 0 ; n1 < node1->nneighbors ; ++n1)
      if( node1->neighbors[n1] == node2 )
         break;
   if( n1 < node1->nneighbors)
      for( n2 = 0 ; n2 < node2->nneighbors ; ++n2)
         if( node2->neighbors[n2] == node1 )
            break;
   if( n1 < node1->nneighbors )
   {
      /* if node2 is neighbor of node1, then node1 is neighbor of node2 */
      assert(node1->neighbors[n1] == node2);
      assert(n2 < node2->nneighbors);
      assert(node2->neighbors[n2] == node1);
      assert(node1->edgeweights[n1] == node2->edgeweights[n2]); /*lint !e777*/
   }
  
   if( n1 == node1->nneighbors || SCIPisLT(scip, weight, node1->edgeweights[n1]) )
   {
      node1->neighbors[n1] = node2;
      node1->edgeweights[n1] = weight;
      node1->relatedrows[n1] = relatedrow;
      node1->nneighbors++;

      node2->neighbors[n2] = node1;
      node2->edgeweights[n2] = weight;
      node2->relatedrows[n2] = relatedrow;
      node2->nneighbors++;

      node1copy->neighbors[n1] = node2copy;
      node1copy->edgeweights[n1] = weight;
      node1copy->relatedrows[n1] = relatedrow;
      node1copy->nneighbors++;

      node2copy->neighbors[n2] = node1copy;
      node2copy->edgeweights[n2] = weight;
      node2copy->relatedrows[n2] = relatedrow;
      node2copy->nneighbors++;    
   }

   return SCIP_OKAY;
}


/** Dijkstra's shortest path algorithm. Calculates the shortest path between
    sourcenode and targetnode. The calculation is aborted if the shortest path
    cannot be shorter than maxdistance */ 
static
SCIP_RETCODE dijkstra(
   SCIP*                 scip,               /**< SCIP data structure */
   ZEROHALF_AUXGRAPH*    graph,              /**< auxiliary graph */
   ZEROHALF_AUXGRAPH_NODE* sourcenode,       /**< start node */
   ZEROHALF_AUXGRAPH_NODE* targetnode,       /**< end node */
   SCIP_Real             maxdistance         /**< calculation will be aborted if a proof is found that no shortest path with
                                              *   length less than maxdistance exists */
   )
{
   ZEROHALF_AUXGRAPH_NODE**   unprocessednodes;
   int                        nunprocessednodes;  
   int                        u;
   int                        v;
   int                        n;
   SCIP_Real                  mindistance;  
   SCIP_Real                  newdistance;
   ZEROHALF_AUXGRAPH_NODE*    currentnode;
  
  
   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->nnodes > 0);
   assert(sourcenode != NULL);
   assert(targetnode != NULL);
   assert(maxdistance > 0.0);
   assert(maxdistance <= 1.0);
  
   /* allocate temporary memory */
   SCIP_CALL(SCIPallocBufferArray(scip, &unprocessednodes, 2 * graph->nnodes));

   /* initialize */
   nunprocessednodes = 0;
   mindistance = 1.0;  
   for( v = 0; v < graph->nnodes ; ++v)
   { 
      graph->nodes[v]->distance              =  1.0;
      graph->nodes[v]->previous              = NULL;
    
      graph->nodecopies[v]->distance         =  1.0;    
      graph->nodecopies[v]->previous         = NULL;

      unprocessednodes[nunprocessednodes]    = graph->nodes[v];
      ++nunprocessednodes;
      unprocessednodes[nunprocessednodes]    = graph->nodecopies[v];
      ++nunprocessednodes;
   }  
   sourcenode->distance = 0.0;
   sourcenode->previous = NULL;
  
   assert(nunprocessednodes == 2 * graph->nnodes);
   assert(nunprocessednodes > 0);

   /* for all nodes */
   while( nunprocessednodes > 0 )
   {
      /* get unprocessed node with minimum distance from sourcenode */
      u = 0;
      mindistance = unprocessednodes[0]->distance;
      for( v = 1 ; v < nunprocessednodes ; ++v)
         if( unprocessednodes[v]->distance < mindistance )
         {
            u = v;
            mindistance = unprocessednodes[v]->distance;
         }
      /* if mindistance is greater than maxdistance then abort dijkstra */
      if( SCIPisGT(scip, mindistance, maxdistance) )
         goto exitdijkstra;
      /* set minimum distance node as currentnode */
      currentnode = unprocessednodes[u];
      unprocessednodes[u] = unprocessednodes[nunprocessednodes-1];
      nunprocessednodes--;

      /* for all neighbors of currentnode */
      for( n = 0 ; n < currentnode->nneighbors ; ++n)
      {
         newdistance = currentnode->distance + currentnode->edgeweights[n];
         /* check if distance via currentnode is less than their actual distance */
         if( SCIPisLT(scip, newdistance , currentnode->neighbors[n]->distance) )
         {
            currentnode->neighbors[n]->distance = newdistance;
            currentnode->neighbors[n]->previous = currentnode;
            /* if neighbors is targetnode then abort dijkstra */
            if( currentnode->neighbors[n] == targetnode )
               goto exitdijkstra;
         }
      }    
   }
  
 exitdijkstra:
   SCIPfreeBufferArray(scip, &unprocessednodes);
   return SCIP_OKAY;
}


/**  separates violated zerohalf cuts by searching for minweight odd-valued
     cycles within an auxiliary graph. (exact method, but only applicable if
     each row contains at most two odd entries; polynomial time) */
static
SCIP_RETCODE separateByAuxGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */        
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered (preprocessed) subproblem mod 2 */
   char                  normtype,           /**< SCIP normtype */
   int                   maxsepacuts,        /**< maximal number of zerohalf cuts separated per separation round */
   int                   maxcuts,            /**< maximal number of zerohalf cuts found per separation round (incl. ineff. cuts) */
   int*                  nsepacuts,          /**< pointer to store current number of separated zerohalf cuts */
   int*                  nzerohalfcuts,      /**< pointer to store current number of found zerohalf cuts */ 
   ZEROHALF_CUTDATA**    zerohalfcuts,       /**< array to store found zerohalf cuts */
   SCIP_Real**           varsolvals,         /**< dense variable LP solution vector */
   SCIP_RESULT*          result,             /**< pointer to SCIP result value of separation */
   SCIP_Bool*            wrongstructure      /**< pointer to store if there is a row with more than two odd entries */
   )
{
   ZEROHALF_AUXGRAPH*    auxgraph;
   int                   rowsindex;
   int                   i;
   int                   j;  
   int                   k;
   int                   n;
   int                   q;

   SCIP_Real*            weights;
   int                   nrowsincut;
   int                   nrrowsincut;
   BITARRAY              rrowsincut;
   ZEROHALF_AUXGRAPH_NODE* node;
  
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(lpdata != NULL);
   assert(mod2data != NULL);
   assert(maxsepacuts >= 0);
   assert(maxcuts >= 0);  
   assert(nsepacuts != NULL);
   assert(nzerohalfcuts != NULL);
   assert(zerohalfcuts != NULL);
   assert(*nsepacuts <= *nzerohalfcuts);
   assert(varsolvals != NULL);
   assert(result != NULL);
   assert(wrongstructure != NULL);

   *wrongstructure = FALSE;
  
   /* check if( A mod 2, b mod 2) is empty */
   if( mod2data->nrows == 0 || mod2data->nrowsind == 0 )
      return SCIP_OKAY;

   /* check if enough cuts have been found */
   if( *nsepacuts >= maxsepacuts || *nzerohalfcuts >= maxcuts )
      return SCIP_OKAY; 

   /* check if matrix has proper structure */ 
   if( !hasMatrixMax2EntriesPerRow(mod2data) )
   {
      *wrongstructure = TRUE;
      return SCIP_OKAY;
   }
  
   /* check if only one row exists */
   if( mod2data->nrowsind == 1 )
   {
      SCIP_CALL(separateByEnumerationHeuristics(scip, sepa, sepadata, lpdata, mod2data,
            normtype, maxsepacuts, maxcuts, nsepacuts, nzerohalfcuts, zerohalfcuts,
            varsolvals, result, 1));
      return SCIP_OKAY;
   }

  
   /* build auxiliary graph */
   SCIP_CALL(ZerohalfAuxGraphCreate(scip, &auxgraph));

   /* create nodes */
   auxgraph->nnodes = mod2data->ncolsind + 1;  
   SCIP_CALL(SCIPallocMemoryArray(scip, &(auxgraph->nodes), auxgraph->nnodes));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(auxgraph->nodecopies), auxgraph->nnodes));
   q = auxgraph->nnodes - 1;
  
   for( j = 0 ; j < auxgraph->nnodes ; ++j)
   {
      SCIP_CALL(ZerohalfAuxGraphNodeCreate(scip, &(auxgraph->nodes[j])));
      SCIP_CALL(ZerohalfAuxGraphNodeCreate(scip, &(auxgraph->nodecopies[j])));
   }
   /* create edges */
   for( i = 0 ; i < mod2data->nrowsind ; ++i)
   {
      rowsindex = mod2data->rowsind[i];
      for( j = 0 ; j < mod2data->ncolsind; ++j)
         if( BITARRAYBITISSET(mod2data->rows[rowsindex], mod2data->colsind[j]) ) /*lint !e701*/
            break;
      for( k = j+1 ; k < mod2data->ncolsind; ++k)
         if( BITARRAYBITISSET(mod2data->rows[rowsindex], mod2data->colsind[k]) ) /*lint !e701*/
            break;

      /* check if row i is a zero row */
      if( j >= mod2data->ncolsind )
      {
         if( mod2data->rhs[rowsindex] )
            if( SCIPisLE(scip, mod2data->slacks[rowsindex], sepadata->maxslack) )
            {
               /* violated {0,1/2} cut has been found */
               weights = NULL;
               nrrowsincut = 1;
               nrowsincut = 0;
               SCIP_CALL(getZerohalfWeightvectorFromSelectedRowsBitarray(scip, sepadata, lpdata, mod2data,
                     mod2data->rowaggregations[rowsindex], &weights, &nrowsincut));
               if( weights == NULL )
               {
                  continue;
               }
               assert(nrowsincut > 0);

               /* create zerohalf cut */         
               SCIP_CALL(ZerohalfCutDataCreate(scip, &(zerohalfcuts[*nzerohalfcuts]),
                     mod2data->relatedsubproblem, mod2data, nrrowsincut, nrowsincut, AUXGRAPH));
               SCIP_CALL(createZerohalfCutFromZerohalfWeightvector(scip, sepa, sepadata,
                     lpdata, weights, normtype, *nzerohalfcuts, varsolvals, zerohalfcuts[*nzerohalfcuts]));
                    
               /* add cut */
               SCIP_CALL( addZerohalfCutToLP(scip, sepadata, zerohalfcuts[*nzerohalfcuts], nsepacuts, result) );
               (*nzerohalfcuts)++;
          
               /* free temporary memory */
               assert( weights != NULL );
               SCIPfreeMemoryArray(scip, &weights);
               weights = NULL;
               
            }
         continue;
      }
    
      /* check if row i has only one entry */
      if( k >= mod2data->ncolsind )
      {
         /* add edges (j,q), (j',q')  or  (j,q'), (j',q)  w.r.t. rhs_i mod 2 */
         SCIP_CALL(addEdgeToAuxGraph(scip, auxgraph, j, q,
               mod2data->rhs[rowsindex],  mod2data->slacks[rowsindex], rowsindex));
         continue;
      }
    
      /* row i has two entries */
      /*   add edges (j,k), (j',k')  or  (j,k'), (j',k)  w.r.t. rhs_i mod 2 */
      SCIP_CALL(addEdgeToAuxGraph(scip, auxgraph, j, k,
            mod2data->rhs[rowsindex], mod2data->slacks[rowsindex], rowsindex));
      
   }
   /* create edges (j,q) and (j',q') for all nodes j */
   for( n = 0 ; n < auxgraph->nnodes && n != q ; ++n)
   {
      SCIP_CALL(addEdgeToAuxGraph(scip, auxgraph, n, q,
            FALSE, mod2data->fracsol[mod2data->colsind[n]], -1));
   }
  
   if( auxgraph->nnodes == 0 )
   {
      /* free temporary memory */
      SCIP_CALL(ZerohalfAuxGraphFree(scip, &auxgraph));
      return SCIP_OKAY;
   }

   weights = NULL;
  
   /* calculate shortest (node_i, nodecopy_i)-paths using the dijkstra algorithm */
   for( n = 0 ; n < auxgraph->nnodes && n != q ; ++n)
   {
      SCIP_CALL(dijkstra(scip, auxgraph,
            auxgraph->nodes[n], auxgraph->nodecopies[n], sepadata->maxslack));

      if( SCIPisLE(scip, auxgraph->nodecopies[n]->distance, sepadata->maxslack) )
      {
         /* a violated {0,1/2} cut has been found */

         /* determine original rows that have to be combined */
         SCIP_CALL(SCIPallocBufferArray(scip, &rrowsincut, mod2data->rowaggregationsbitarraysize));
         BITARRAYCLEAR(rrowsincut, mod2data->rowaggregationsbitarraysize);

         nrrowsincut = 0;
         nrowsincut = 0;
         node = auxgraph->nodecopies[n];
         while( node->previous != NULL )
         {
            for( i = 0 ; i < node->nneighbors ; ++i)
               if( node->neighbors[i] == node->previous )
                  break;
            assert(i < node->nneighbors);
            rowsindex = node->relatedrows[i];
            BITARRAYSXOR(mod2data->rowaggregations[rowsindex], rrowsincut ,
               mod2data->rowaggregationsbitarraysize);
            node = node->previous;
            ++nrrowsincut;
         }

         /* get {0,1/2}-weight vector */
         assert(weights == NULL);
         SCIP_CALL(getZerohalfWeightvectorFromSelectedRowsBitarray(scip, sepadata, lpdata,
               mod2data, rrowsincut, &weights, &nrowsincut));
         if( weights == NULL )
         { 
            SCIPfreeBufferArray(scip, &rrowsincut);
            continue;
         }
         assert(nrowsincut > 0);
      
         /* create zerohalf cut */
         SCIP_CALL(ZerohalfCutDataCreate(scip, &(zerohalfcuts[*nzerohalfcuts]),
               mod2data->relatedsubproblem, mod2data, nrrowsincut, nrowsincut, AUXGRAPH));
         SCIP_CALL(createZerohalfCutFromZerohalfWeightvector(scip, sepa, sepadata,
               lpdata, weights, normtype, *nzerohalfcuts, varsolvals, zerohalfcuts[*nzerohalfcuts]));

         /* add cut */
         SCIP_CALL( addZerohalfCutToLP(scip, sepadata, zerohalfcuts[*nzerohalfcuts], nsepacuts, result) );
         (*nzerohalfcuts)++;

         /* free temporary memory */
         assert( weights != NULL );
         SCIPfreeMemoryArray(scip, &weights);
         weights = NULL;
         
         if( rrowsincut != NULL )
         {
            SCIPfreeBufferArray(scip, &rrowsincut);
            rrowsincut = NULL;
         }
      }
   }
  
   /* free temporary memory */
   SCIP_CALL(ZerohalfAuxGraphFree(scip, &auxgraph));

   return SCIP_OKAY;
}


 
/** separates violated zerohalf cuts using an extended Gaussian elimination. (heuristic; polynomial time) */
static
SCIP_RETCODE separateByGaussHeuristics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */        
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered (preprocessed) subproblem mod 2 */
   char                  normtype,           /**< SCIP normtype */
   int                   maxsepacuts,        /**< maximal number of zerohalf cuts separated per separation round */
   int                   maxcuts,            /**< maximal number of zerohalf cuts found per separation round (incl. ineff. cuts) */
   int*                  nsepacuts,          /**< pointer to store current number of separated zerohalf cuts */
   int*                  nzerohalfcuts,      /**< pointer to store current number of found zerohalf cuts */
   ZEROHALF_CUTDATA**    zerohalfcuts,       /**< array to store found zerohalf cuts */
   SCIP_Real**           varsolvals,         /**< dense variable LP solution vector */
   SCIP_RESULT*          result              /**< pointer to SCIP result value of separation */
   )
{
   int                   pivotrow;
   int                   pivotcol;
   int                   identsubmatrixsize;
   int                   rowsbind;
   BITARRAYBITMASKTYPE   rowsbmask;
   int                   r;
   int                   temp;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(lpdata != NULL);
   assert(mod2data != NULL);
   assert(maxsepacuts >= 0);
   assert(maxcuts >= 0);  
   assert(nsepacuts != NULL);
   assert(nzerohalfcuts != NULL);
   assert(zerohalfcuts != NULL);
   assert(*nsepacuts <= *nzerohalfcuts);
   assert(varsolvals != NULL);
   assert(result != NULL);

  
   /* check if( A mod 2, b mod 2) is empty */
   if( mod2data->nrows == 0 || mod2data->nrowsind == 0 )
      return SCIP_OKAY;

   /* check if enough cuts have been found */
   if( *nsepacuts >= maxsepacuts || *nzerohalfcuts >= maxcuts )
      return SCIP_OKAY; 

   if( mod2data->nrowsind == 1 )
   {
      SCIP_CALL(separateByEnumerationHeuristics(scip, sepa, sepadata, lpdata, mod2data,
            normtype, maxsepacuts, maxcuts, nsepacuts, nzerohalfcuts, zerohalfcuts,
            varsolvals, result, 1));
      return SCIP_OKAY;
   }

   identsubmatrixsize = 0;
    
   /* apply Gaussian elimination mod 2 */

   /* choose pivot col */
   for( pivotcol = 0; pivotcol < mod2data->ncolsind; ++pivotcol )
   {
      if( identsubmatrixsize == mod2data->nrowsind )
         break;
    
      /* sort row indices sets w.r.t. to their slack values NON-DECREASINGLY */
      SCIPsortInd(mod2data->rowsind + identsubmatrixsize , compRealNonDecreasing , 
         (void*) mod2data->slacks , mod2data->nrowsind - identsubmatrixsize);
    
      /* break if no unprocessed row with slack <= maxslack is left */
      if( SCIPisGT(scip, mod2data->slacks[mod2data->rowsind[identsubmatrixsize]], sepadata->maxslack) )
         break;
    
      /* determine pivot row */
      rowsbind = (int) GETBITARRAYINDEX(mod2data->colsind[pivotcol]);
      rowsbmask = GETBITARRAYMASK(mod2data->colsind[pivotcol]); /*lint !e701*/
      for( pivotrow = identsubmatrixsize ; pivotrow < mod2data->nrowsind ; ++pivotrow)
         if( mod2data->rows[mod2data->rowsind[pivotrow]][rowsbind] & rowsbmask )
            break;
      if( pivotrow == mod2data->nrowsind )
         continue;

      /* Gaussian elimination step */
      for( r = 0 ; r < mod2data->nrowsind
              && SCIPisLE(scip, mod2data->slacks[mod2data->rowsind[r]], sepadata->maxslack) ; ++r) 
      {
         if( r == pivotrow )
            continue;
         if( mod2data->rows[mod2data->rowsind[r]][rowsbind] & rowsbmask )
         {
            /* add pivot row to r-th row */
            mod2data->slacks[mod2data->rowsind[r]] += mod2data->slacks[mod2data->rowsind[pivotrow]];

#ifndef SCIP_DEBUG
            /* avoid expensive operations on rows with slack > maxslack */
            if( SCIPisLE(scip, mod2data->slacks[mod2data->rowsind[r]], sepadata->maxslack) )
#endif
            {
               BITARRAYSXOR(mod2data->rows[mod2data->rowsind[pivotrow]],
                  mod2data->rows[mod2data->rowsind[r]],mod2data->rowsbitarraysize);
               BITARRAYSXOR(mod2data->rowaggregations[mod2data->rowsind[pivotrow]],
                  mod2data->rowaggregations[mod2data->rowsind[r]],mod2data->rowaggregationsbitarraysize);
               mod2data->rhs[mod2data->rowsind[r]] =
                  XOR(mod2data->rhs[mod2data->rowsind[pivotrow]],mod2data->rhs[mod2data->rowsind[r]]);
            }
         }      
      }
    
      /* swap index set positions */
      temp = mod2data->rowsind[pivotrow];
      mod2data->rowsind[pivotrow] = mod2data->rowsind[identsubmatrixsize];
      mod2data->rowsind[identsubmatrixsize] = temp;
      temp = mod2data->colsind[pivotcol];
      mod2data->colsind[pivotcol] = mod2data->colsind[identsubmatrixsize];
      mod2data->colsind[identsubmatrixsize] = temp;
    
      identsubmatrixsize++;
   }

   /* remove (generated) column singletons */
   SCIP_CALL(preprocessColumns(scip, sepadata, lpdata, mod2data,
         0, mod2data->ncolsind, FALSE, TRUE, TRUE));

   /* remove zero rows and rows with slack > maxslack */
   SCIP_CALL(preprocessRows(scip, sepadata, lpdata, mod2data,
         0, mod2data->nrowsind, TRUE, TRUE, FALSE));

   /* search for zerohalf cuts */ 
   SCIP_CALL(preprocessTrivialZerohalfCuts(scip, sepa, sepadata, lpdata, mod2data,
         0, mod2data->nrowsind, normtype, maxsepacuts, maxcuts, nsepacuts,
         nzerohalfcuts, zerohalfcuts, varsolvals, HEURISTICSGAUSS, result));
  
   return SCIP_OKAY;  
}


/** processes subproblem (i.e. runs separation algorithms)*/
static
SCIP_RETCODE process(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */    
   ZEROHALF_LPDATA*      lpdata,             /**< data of current LP relaxation */
   ZEROHALF_MOD2DATA*    mod2data,           /**< considered (preprocessed) subproblem mod 2 */
   char                  normtype,           /**< SCIP normtype */
   int                   maxsepacuts,        /**< maximal number of zerohalf cuts separated per separation round */
   int                   maxcuts,            /**< maximal number of zerohalf cuts found per separation round (incl. ineff. cuts) */
   int*                  nsepacuts,          /**< pointer to store current number of separated zerohalf cuts */
   int*                  nzerohalfcuts,      /**< pointer to store current number of found zerohalf cuts */
   ZEROHALF_CUTDATA**    zerohalfcuts,       /**< array to store found zerohalf cuts */
   SCIP_Real**           varsolvals,         /**< dense variable LP solution vector */
   SCIP_RESULT*          result              /**< pointer to SCIP result value of separation */
   )
{  
   char                  sepamethod;
   int                   i;
   SCIP_Bool             stop;

   char                  sepaname[SCIP_MAXSTRLEN];
#ifdef ZEROHALF__PRINT_STATISTICS
   int                   nsepacutsbefore;
   int                   nsepacutsinitial;
   int                   nzerohalfcutsbefore;
   int                   nzerohalfcutsinitial;
   SCIP_CLOCK* timer;
   SCIP_CLOCK* sepatimer;
#endif

   int                   ncutsfoundbefore;
   SCIP_Bool             wrongstructure;
  
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(lpdata != NULL);
   assert(mod2data != NULL);
   assert(maxsepacuts >= 0);
   assert(maxcuts >= 0);  
   assert(result != NULL);
   assert(nsepacuts != NULL);
   assert(nzerohalfcuts != NULL);
   assert(zerohalfcuts != NULL);
   assert(*nsepacuts <= *nzerohalfcuts);
  
   assert(mod2data->relatedsubproblem != NULL);
   assert(mod2data->rows != NULL);
   assert(mod2data->rowaggregations != NULL);
   assert(mod2data->rhs != NULL);
   assert(mod2data->slacks != NULL);
   assert(mod2data->fracsol != NULL);
   assert(mod2data->nrows > 0);
   assert(mod2data->rowsind != NULL);
   assert(mod2data->colsind != NULL);

#ifdef ZEROHALF__PRINT_STATISTICS
   nsepacutsinitial = *nsepacuts;
   nzerohalfcutsinitial = *nzerohalfcuts;
#endif

   if( sepadata->nsepamethods == -1 )
   {
      sepadata->nsepamethods = (int) strlen(sepadata->sepamethods);
      if( sepadata->nsepamethods > 0 && sepadata->sepamethods[0] == '-' )
         sepadata->nsepamethods = 0;    
   }
  
   if( sepadata->nsepamethods == 0 )
      return SCIP_OKAY;

  
   /* statistics */  
#ifdef ZEROHALF__PRINT_STATISTICS
   if( sepadata->sepatimers == NULL )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(sepadata->sepatimers), sepadata->nsepamethods + 1) );
      for( i = 0 ; i < sepadata->nsepamethods + 1 ; ++i)
      {
         ZEROHALFcreateTimer((sepadata->sepatimers[i]));
      }
      SCIP_CALL( SCIPallocMemoryArray(scip, &(sepadata->nsepacutsalgo), sepadata->nsepamethods + 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(sepadata->nzerohalfcutsalgo), sepadata->nsepamethods + 1) );
      BMSclearMemoryArray(sepadata->nsepacutsalgo, sepadata->nsepamethods + 1);
      BMSclearMemoryArray(sepadata->nzerohalfcutsalgo, sepadata->nsepamethods + 1);      
   }
 
   nsepacutsinitial = *nsepacuts;
   nzerohalfcutsinitial = *nzerohalfcuts;
    
   ZEROHALFstatisticsMessage("\n");
   ZEROHALFstatisticsMessage("                | -------------------------- subproblem ----\
---------- | - callback (algo) - | ----- callback ---- | --total-\n");
   ZEROHALFstatisticsMessage("                | nrowsind | ncolsind | nsepcuts | ncutsfnd \
|     time | nsepcuts | ncutsfnd | nsepcuts | ncutsfnd |     time\n");
   ZEROHALFstatisticsMessage("%15s | %8d | %8d | %8d | %8d | %8.4f | %8s | %8s | %8d | %8d | %8.4f\n",
      "START PROCESSG", mod2data->nrowsind, mod2data->ncolsind,
      0, 0, 0.0, " ", " ", *nsepacuts, *nzerohalfcuts, ZEROHALFevalTimer(sepadata->sepatimers[sepadata->nsepamethods]));
#endif

   if( mod2data->nrowsind == 0 || mod2data->ncolsind == 0 )
      return SCIP_OKAY;

#ifdef ZEROHALF__PRINT_STATISTICS
   ZEROHALFcreateNewTimer(timer);
   ZEROHALFcreateNewTimer(sepatimer);  
   ZEROHALFstartTimer(sepadata->sepatimers[sepadata->nsepamethods]);  
   ZEROHALFstartTimer(timer);
#endif

   stop = FALSE;
   for( i = 0 ; i < sepadata->nsepamethods && !stop; ++i)
   {
      /* statistics*/ 
#ifdef ZEROHALF__PRINT_STATISTICS
      ZEROHALFstartTimer(sepatimer); 
      ZEROHALFstartTimer(sepadata->sepatimers[i]);
      nsepacutsbefore = *nsepacuts;
      nzerohalfcutsbefore = *nzerohalfcuts;
    
      /* abort if enough cuts have already been found */
      if( *nsepacuts >= maxsepacuts || *nzerohalfcuts >= maxcuts )
         break;
      /* abort if matrix is empty */
      if( mod2data->nrowsind == 0 || mod2data->ncolsind == 0 )
         break;
#endif

      /* process */
      sepamethod = sepadata->sepamethods[i];    

      switch(sepamethod)
      {
      case STOPIFCUTWASFOUND:
         if( *nsepacuts > 0 )
         {
            stop = TRUE;
         }      
         strncpy(sepaname,"nsepcuts>0 ?",SCIP_MAXSTRLEN);
         break;
      case SOLVEAUXSCIP: 
      case SOLVEAUXSCIPEXACT:
         ncutsfoundbefore = *nzerohalfcuts;
         SCIP_CALL(separateBySolvingAuxIP(scip, sepa, sepadata, lpdata, mod2data, normtype,
               maxsepacuts, maxcuts, (sepamethod == SOLVEAUXSCIP),
               nsepacuts, nzerohalfcuts, zerohalfcuts, varsolvals, result));/*lint !e641*/
         strncpy(sepaname,"auxiliary ip",SCIP_MAXSTRLEN);
         if( *nzerohalfcuts == ncutsfoundbefore )
         {
            /* no violated cut has been found. hence a proof of non-existence is given */
            stop = TRUE;
         }        
         break;
      case ENUMHEURNMAX1:
         SCIP_CALL(separateByEnumerationHeuristics(scip, sepa, sepadata, lpdata, mod2data, normtype, maxsepacuts,
               maxcuts, nsepacuts, nzerohalfcuts, zerohalfcuts, varsolvals, result, 1));
         strncpy(sepaname,"enum k=1",SCIP_MAXSTRLEN);
         break;
      case ENUMHEURNMAX2:
         SCIP_CALL(separateByEnumerationHeuristics(scip, sepa, sepadata, lpdata, mod2data, normtype, maxsepacuts,
               maxcuts, nsepacuts, nzerohalfcuts, zerohalfcuts, varsolvals, result, 2));
         strncpy(sepaname,"enum k=1..2",SCIP_MAXSTRLEN);
         break;
      case GAUSSHEUR:
         SCIP_CALL(separateByGaussHeuristics(scip, sepa, sepadata, lpdata, mod2data, normtype, maxsepacuts, maxcuts,
               nsepacuts, nzerohalfcuts, zerohalfcuts, varsolvals, result));
         strncpy(sepaname,"Gauss heur",SCIP_MAXSTRLEN);
         break;
      case MAX2ODDENTRIESPERROW:
         ncutsfoundbefore = *nzerohalfcuts;
         SCIP_CALL(separateByAuxGraph(scip, sepa, sepadata, lpdata, mod2data, normtype, maxsepacuts, maxcuts,
               nsepacuts, nzerohalfcuts, zerohalfcuts, varsolvals, result, &wrongstructure));
         strncpy(sepaname,"auxgraph",SCIP_MAXSTRLEN);
         if( ! wrongstructure )
         {
            if( *nzerohalfcuts == ncutsfoundbefore )
            {
               /* no violated cut has been found. hence a proof of non-existence is given */
               stop = TRUE;
            }
         }
         break;
      default:
         SCIPerrorMessage("invalid sepamethod '%c'\n", sepadata->sepamethods[i]);
         return SCIP_INVALIDDATA;   
      }

#ifdef ZEROHALF__PRINT_STATISTICS
      /* statistics */
      ZEROHALFstopTimer(sepadata->sepatimers[i]);
      ZEROHALFstopTimer(sepatimer);
      sepadata->nsepacutsalgo[i] += *nsepacuts - nsepacutsbefore;
      sepadata->nzerohalfcutsalgo[i] += *nzerohalfcuts - nzerohalfcutsbefore;
      ZEROHALFstatisticsMessage("%15s | %8d | %8d | %8d | %8d | %8.4f | %8d | %8d | %8d | %8d | %8.4f\n",
         sepaname, mod2data->nrowsind, mod2data->ncolsind,
         *nsepacuts - nsepacutsbefore, *nzerohalfcuts - nzerohalfcutsbefore,
         ZEROHALFevalTimer(sepatimer),
         sepadata->nsepacutsalgo[i], sepadata->nzerohalfcutsalgo[i],
         *nsepacuts, *nzerohalfcuts, ZEROHALFevalTimer(sepadata->sepatimers[i]));
      ZEROHALFresetTimer(sepatimer); 
#endif
   }

#ifdef ZEROHALF__PRINT_STATISTICS 
   /* statistics */
   ZEROHALFstopTimer(timer);
   ZEROHALFstatisticsMessage("%15s | %8d | %8d | %8d | %8d | %8.4f | %8s | %8s | %8d | %8d | %8.4f\n",
      "PROCESSED", mod2data->nrowsind, mod2data->ncolsind,
      *nsepacuts - nsepacutsinitial , *nzerohalfcuts - nzerohalfcutsinitial,
      ZEROHALFevalTimer(timer), " ", " ", *nsepacuts, *nzerohalfcuts,
      ZEROHALFevalTimer(sepadata->sepatimers[sepadata->nsepamethods])); 
   ZEROHALFfreeTimer(timer);
   ZEROHALFfreeTimer(sepatimer);
   ZEROHALFstopTimer(sepadata->sepatimers[sepadata->nsepamethods]);
#endif

   return SCIP_OKAY;
}


#ifdef ZEROHALF__PRINT_STATISTICS 
/** prints statistical information about the found zerohalfcuts as table */ 
static
void printZerohalfCutsStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */     
   ZEROHALF_CUTDATA**    zerohalfcuts,       /**< array of zerohalf cuts */
   int                   nzerohalfcuts,      /**< number of zerohalf cuts */
   int*                  zerohalfcutsindices,/**< sorted index set (or NULL) */
   SCIP_Real*            zerohalfcutsprios,  /**< sorted priorities (or NULL) */
   SCIP_Real*            zerohalfcutsminortho,/**< sorted minimal orthogonalities (or NULL) */
   int                   nsepacuts           /**< number of separated zerohalf cuts */
   )
{
   ZEROHALF_CUTDATA*     cut;
   int                   si;
   int                   i;

   assert(scip != NULL);
   assert(zerohalfcuts != NULL);
   assert(nzerohalfcuts >= 0);

   ZEROHALFstatisticsMessage("\n");
   ZEROHALFstatisticsMessage(" NZEROHALFCUTS:   %6d\n", nzerohalfcuts);
   ZEROHALFstatisticsMessage(" NSEPACUTS:       %6d\n", nsepacuts);

   if( nzerohalfcuts == 0 )
      return;
  
   ZEROHALFstatisticsMessage("%15s |  index | A |     viol | efficacy | ef? \
| minortho |    #nonz |     norm | #origrows | #preprows | by | priority\n", " ");
   for( i = 0 ; i < nzerohalfcuts ; ++i)
   {  
      if( zerohalfcutsindices == NULL )
         si = i;
      else
         si = zerohalfcutsindices[i];
      assert(0 <= i);
      assert(i < nzerohalfcuts);
      cut = zerohalfcuts[si];
      ZEROHALFstatisticsMessage("%15s | %6d | %1s | %8f | %8f | %3s | %8f | %8d | %8.2f | %9d | %9d | %2s | %8.6f\n", 
         "ZEROHALFCUT", si,
         cut->addedtolp ? "A" : "-",
         cut->violation, cut->efficacy,
         SCIPisEfficacious(scip, cut->efficacy) ? "yes" : " no",
         zerohalfcutsminortho != NULL ? zerohalfcutsminortho[i] : 0.0,
         cut->nnonz, cut->norm, cut->nrowsincut, cut->nrrowsincut,
         (cut->separatedby == AUXIP ? "IP" :
            (cut->separatedby == DECOMPOSITION ? "DC" :
               (cut->separatedby == PPZEROONEROW ? "PP" :
                  (cut->separatedby == HEURISTICSENUM ? "HE" :
                     (cut->separatedby == HEURISTICSGAUSS ? "HG" :
                        (cut->separatedby == AUXGRAPH ? "2R" :
                           "?"
                           )))))),
         zerohalfcutsprios != NULL ? zerohalfcutsprios[i] : 0.0
         );      
   }
}
#endif


 
/* --------------------------------------------------------------------------------------------------------------------
 * callback methods of separator
 * -------------------------------------------------------------------------------------------------------------------- */


/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyZerohalf)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaZerohalf(scip) );
 
   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeZerohalf)
{  
   SCIP_SEPADATA* sepadata;

   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
  
   if( sepadata->pptimers != NULL )
   {
#ifdef ZEROHALF__PRINT_STATISTICS   
      int i;

      for( i = 0 ; i < sepadata->nppmethods + 1 ; ++i)
      {
         ZEROHALFfreeTimer((sepadata->pptimers[i]));
      }
#endif
      SCIPfreeMemoryArray(scip, &(sepadata->pptimers));
   }

   if( sepadata->sepatimers != NULL )
   {
#ifdef ZEROHALF__PRINT_STATISTICS   
      int i;

      for( i = 0 ; i < sepadata->nsepamethods + 1 ; ++i)
      {
         ZEROHALFfreeTimer((sepadata->sepatimers[i]));
      }
#endif
      SCIPfreeMemoryArray(scip, &(sepadata->sepatimers));
   }

#ifdef ZEROHALF__PRINT_STATISTICS   
   if( sepadata->dtimer != NULL )
   {
      ZEROHALFfreeTimer((sepadata->dtimer));
   }  
#endif

   if( sepadata->nsepacutsalgo != NULL )
   {
      SCIPfreeMemoryArray(scip, &(sepadata->nsepacutsalgo));
   }
   if( sepadata->nzerohalfcutsalgo != NULL )
   {
      SCIPfreeMemoryArray(scip, &(sepadata->nzerohalfcutsalgo));
   }
   if( sepadata->origrows != NULL )
   {
      SCIPfreeMemoryArray(scip, &(sepadata->origrows));
   }
   SCIPfreeMemory(scip, &sepadata);   
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpZerohalf)
{ 
   SCIP_SEPADATA*        sepadata;
   ZEROHALF_LPDATA*      lpdata;
   ZEROHALF_MOD2DATA*    mod2data;
   ZEROHALF_CUTDATA**    zerohalfcuts;  
   SCIP_Real*            varsolvals;    
   SCIP_Real*            subproblempriorities;  
   int*                  sortedsubproblems;  
   int                   ncalls;
   int                   depth;
   int                   maxsepacuts;
   int                   maxcuts;
   int                   nsepacuts;
   int                   nzerohalfcuts;
   int                   subproblemindex;
   int                   i;
   int                   j;
   char                  normtype;
   SCIP_Longint          totalncalls;
   SCIP_Real             nrcolsfactor;

   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;  
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   depth = SCIPgetDepth(scip); 
   ncalls = SCIPsepaGetNCallsAtNode(sepa); 
   totalncalls = SCIPsepaGetNCalls(sepa);

   /* only call separator, if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* only call the {0,1/2}-cut separator a given number of times at all */
   if( sepadata->maxncalls > -1 && totalncalls > sepadata->maxncalls - 1 )
      return SCIP_OKAY;

   /* only call the {0,1/2}-cut separator a given number of times at each node */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* only call separator if depth<=maxdepth or maxdepth unlimited */
   if( sepadata->maxdepth > -1 && depth > sepadata->maxdepth )
      return SCIP_OKAY;
  
   /* only call separator, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call separator, if there are fractional variables */
   if( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   /* allocate temporary memory for LP data structures */
   SCIP_CALL(ZerohalfLPDataCreate(scip, &lpdata));
  
   /* get variables data */
   SCIP_CALL(SCIPgetVarsData(scip, &(lpdata->vars), &(lpdata->nvars), NULL, NULL, NULL, NULL));
  
   /* get LP data */
   SCIP_CALL(SCIPgetLPColsData(scip, &(lpdata->cols), &(lpdata->ncols)));
   SCIP_CALL(SCIPgetLPRowsData(scip, &(lpdata->rows), &(lpdata->nrows)));
   if( lpdata->ncols == 0 || lpdata->nrows == 0 )
   {
      SCIP_CALL(ZerohalfLPDataFree(scip, &lpdata));
      return SCIP_OKAY;
   }

   /* store original LP rows indices if necessary */
   if( sepadata->onlyorigrows && sepadata->origrows == NULL )
   {
      SCIP_Bool issorted;                                                 
      int temp;                                                           

      SCIP_CALL( SCIPallocMemoryArray(scip, &(sepadata->origrows), lpdata->nrows) );
      for( i = 0 ; i < lpdata->nrows ; ++i )
         sepadata->origrows[i] = SCIProwGetIndex(lpdata->rows[i]);
      sepadata->norigrows = lpdata->nrows;

      /* bubble sort indices (s.t. binary search can be used to find an index) */
      do                                                                  
      {                                                                   
         issorted = TRUE;
         for( i = 0 ; i < sepadata->norigrows - 1 ; ++i )
         {
            if( sepadata->origrows[i] > sepadata->origrows[i+1] )
            {
               temp = sepadata->origrows[i];
               sepadata->origrows[i] = sepadata->origrows[i+1];
               sepadata->origrows[i+1] = temp;
               issorted = FALSE;
            }
         }
      }
      while( !issorted );
   }

   /* get the maximal number of cuts allowed in a separation round */
   if( depth == 0 )
   {
      maxcuts = sepadata->maxcutsroot;
      maxsepacuts = MIN(maxcuts, sepadata->maxsepacutsroot);
   }
   else
   {
      maxcuts = sepadata->maxcuts;
      maxsepacuts = MIN(maxcuts, sepadata->maxsepacuts);    
   }

   /* only call the {0,1/2}-cut separator if at least one cut is allowed in a separation round */
   if( maxsepacuts == 0 )
   {
      SCIP_CALL(ZerohalfLPDataFree(scip, &lpdata));
      return SCIP_OKAY;
   }

  
#ifdef ZEROHALF__PRINT_STATISTICS  
   ZEROHALFstatisticsMessage("= SEPA_ZEROHALF ================================================================\
=============================================\n"); 
#endif

   /* allocate further temporary memory */
   SCIP_CALL(SCIPallocMemoryArray(scip, &(lpdata->intscalarsleftrow), lpdata->nrows));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(lpdata->intscalarsrightrow), lpdata->nrows));
   
   /* initialize */
   BMSclearMemoryArray(lpdata->intscalarsleftrow, lpdata->nrows);
   BMSclearMemoryArray(lpdata->intscalarsrightrow, lpdata->nrows);

   *result = SCIP_DIDNOTFIND;
   sepadata->maxslack = 1.0 - 2.0 * sepadata->minviolation;
   nsepacuts = 0;
   nzerohalfcuts = 0;
   varsolvals = NULL;
   SCIP_CALL(SCIPgetCharParam(scip, "separating/efficacynorm", &normtype));

#ifdef ZEROHALF__PRINT_STATISTICS
   if( sepadata->nsepamethods > 0 )
   {
      assert(sepadata->nsepacutsalgo != NULL);
      assert(sepadata->nzerohalfcutsalgo != NULL);    
      BMSclearMemoryArray(sepadata->nsepacutsalgo, sepadata->nsepamethods + 1);
      BMSclearMemoryArray(sepadata->nzerohalfcutsalgo, sepadata->nsepamethods + 1);
   }
#endif

   /* determine maximal number of nonzeros allowed in a zerohalf cut */
   sepadata->maxnnonz = 0;
   for( i = 0 ; i < lpdata->nrows ; i++ )
   {
      sepadata->maxnnonz += SCIProwGetNLPNonz(lpdata->rows[i]); 
   }
   sepadata->maxnnonz = (int) floor( 10.0 * sepadata->maxnnonz / (double) lpdata->nrows);
   sepadata->maxnnonz = (int) floor(MIN(0.1 * lpdata->ncols, sepadata->maxnnonz)) + NNONZOFFSET; 

   /* search for relevant columns */
   SCIP_CALL(getRelevantColumns(scip, lpdata));
   if( lpdata->subproblems[0]->nrcols == 0 )
   {
      SCIP_CALL(ZerohalfLPDataFree(scip, &lpdata));
      return SCIP_OKAY;
   }
  
   /* search for relevant rows */
   SCIP_CALL(getRelevantRows(scip, sepadata, lpdata));
   if( lpdata->subproblems[0]->nrrows == 0 )
   {
      SCIP_CALL(ZerohalfLPDataFree(scip, &lpdata));
      return SCIP_OKAY;
   }

#ifdef ZEROHALF__PRINT_STATISTICS
   SCIP_CALL( printPreprocessingStatistics(scip, lpdata) );
#endif
  
   /* try to decompose problem into subproblems (and delete obviously redundant subproblems)*/
   if( sepadata->decomposeproblem )
   {
      SCIP_CALL( decomposeProblem(scip, sepadata, lpdata) ); 
   }

   /* sort subproblems */
   SCIP_CALL(SCIPallocBufferArray(scip, &sortedsubproblems, lpdata->nsubproblems));
   SCIP_CALL(SCIPallocBufferArray(scip, &subproblempriorities, lpdata->nsubproblems));
   nrcolsfactor = 1.0 / (2.0 * lpdata->ncols);
   for( i = 0 ; i < lpdata->nsubproblems; ++i)
   {
      sortedsubproblems[i] = i;
      subproblempriorities[i] = ((SCIP_Real) lpdata->subproblems[i]->nrrows)
         + nrcolsfactor * ((SCIP_Real) lpdata->subproblems[i]->nrcols);
   }
   SCIPsortRealInt(subproblempriorities, sortedsubproblems, lpdata->nsubproblems);

   /* allocate temporary memory for storing separated zerohalf cuts */
   SCIP_CALL(SCIPallocBufferArray(scip, &zerohalfcuts, maxcuts));
 
   /* process each subproblem */
   for( i = 0; i < lpdata->nsubproblems; ++i )
   {
      /* check if enough cuts have been found */
      if( nsepacuts >= maxsepacuts || nzerohalfcuts >= maxcuts )
         break;      

      subproblemindex = sortedsubproblems[i];
      sepadata->maxslack = 1.0 - 2.0 * sepadata->minviolation; /* ignore previous changes */

      /* if the subproblem consists of a single row a^Tx <= b such that a = 0 and b mod 2 = 1,
       * then a violated zerohalf cut can be generated by multiplying this row with 0.5 and rounding down. 
       */
      if( lpdata->subproblems[subproblemindex]->nrrows == 1 && lpdata->subproblems[subproblemindex]->nrcols == 0 )
      {        
         SCIP_Real* weights;

#ifdef SCIP_DEBUG
         SCIP_CALL( debugPrintLPRowsAndCols(scip, lpdata) );
         SCIPdebugMessage("\n");
         debugPrintSubLpData(scip, lpdata, lpdata->subproblems[subproblemindex]);
         SCIPdebugMessage("\n");
#endif
         /* create weightvector */
         weights = NULL;
         SCIP_CALL( getZerohalfWeightvectorForSingleRow(scip, sepadata, lpdata, lpdata->subproblems[subproblemindex]->rrows[0],
               0, &weights) );
         if( weights == NULL )
            continue;

         /* create zerohalf cut */
         SCIP_CALL( ZerohalfCutDataCreate(scip, &(zerohalfcuts[nzerohalfcuts]),
               lpdata->subproblems[subproblemindex], NULL, 1, 1, DECOMPOSITION) );
         SCIP_CALL( createZerohalfCutFromZerohalfWeightvector(scip, sepa, sepadata,
               lpdata, weights, normtype, nzerohalfcuts, &varsolvals, zerohalfcuts[nzerohalfcuts]) );
      
         /* add cut to LP */
         SCIP_CALL( addZerohalfCutToLP(scip, sepadata, zerohalfcuts[nzerohalfcuts], &nsepacuts, result) );
         nzerohalfcuts++;
      
         /* free weightsvector memory */
         SCIPfreeMemoryArray(scip, &weights);

         if( lpdata->rrowsindexofleftrow[lpdata->subproblems[subproblemindex]->rrows[0]] >= 0 )
            lpdata->rrowsindexofleftrow[lpdata->subproblems[subproblemindex]->rrows[0]] =
               DEFINES_VIOLATED_ZEROHALF_CUT;
         else
            if( lpdata->rrowsindexofrightrow[lpdata->subproblems[subproblemindex]->rrows[0]] >= 0 )
               lpdata->rrowsindexofrightrow[lpdata->subproblems[subproblemindex]->rrows[0]] =
                  DEFINES_VIOLATED_ZEROHALF_CUT;
      
         continue;
      }

      /* check if enough cuts have been found */
      if( nsepacuts >= maxsepacuts || nzerohalfcuts >= maxcuts )
         break; 
    
      /* allocate temporary memory for data (mod 2) structures */  
      SCIP_CALL( ZerohalfMod2DataCreate(scip, &mod2data) );
    
      /* store data (mod 2) */    
      SCIP_CALL( storeMod2Data(scip, sepadata, lpdata, subproblemindex, mod2data) );  
     
      /* preprocess subproblem: reduce problem size and/or separate 'easy' zerohalf cuts */
      SCIP_CALL( preprocess(scip, sepa, sepadata, lpdata, mod2data, normtype, maxsepacuts, maxcuts,
            &nsepacuts, &nzerohalfcuts, zerohalfcuts, &varsolvals, result) );

      if( nsepacuts < maxsepacuts || nzerohalfcuts < maxcuts )
      {
         /* process subproblem: separate violated zerohalf cuts */
         SCIP_CALL( process(scip, sepa, sepadata, lpdata, mod2data, normtype, maxsepacuts, maxcuts,
               &nsepacuts, &nzerohalfcuts, zerohalfcuts, &varsolvals, result) );
      }

      /* free temporary memory */
      SCIP_CALL( ZerohalfMod2DataFree(scip, &mod2data) );

      /* update pointer of zerohalfcuts that has become invalid */
      for( j = 0 ; j < nzerohalfcuts ; ++j)
         zerohalfcuts[j]->relatedmod2data = NULL;    
   }
  
#ifdef ZEROHALF__PRINT_STATISTICS
   if( !sepadata->usezhcutpool ) 
      printZerohalfCutsStatistics(scip, sepadata, zerohalfcuts, nzerohalfcuts, NULL, NULL, NULL, nsepacuts);
#endif
  
   if( sepadata->usezhcutpool )
   {
      ZEROHALF_CUTDATA*   cutdatai;
      ZEROHALF_CUTDATA*   cutdataj;
      SCIP_Real*          zerohalfcutpriorities;
      SCIP_Real*          zerohalfcutminortho;
      int*                sortedzerohalfcuts;
      SCIP_Real           violationbucketsize;    
      SCIP_Bool           hasminorthogonality;    
      int                 si;
      int                 sj;
      int                 nignoredcuts;
    
      /* allocate temporary memory */
      SCIP_CALL(SCIPallocBufferArray(scip, &sortedzerohalfcuts, nzerohalfcuts));
      SCIP_CALL(SCIPallocBufferArray(scip, &zerohalfcutpriorities, nzerohalfcuts));
      SCIP_CALL(SCIPallocBufferArray(scip, &zerohalfcutminortho, nzerohalfcuts));

      /* initialize */
      violationbucketsize = 10.0;
      nignoredcuts = 0;

#if 1 /* old cutpool version */ 
      /* sort zerohalf cutpool w.r.t. the violation (primary) and the density (secondary) */    
      for( i = 0 ; i < nzerohalfcuts; ++i)
      {
         sortedzerohalfcuts[i] = i;
         zerohalfcutpriorities[i] = SCIPfloor(scip, violationbucketsize * (zerohalfcuts[i])->violation)
            + (1.0 - (SCIP_Real) (zerohalfcuts[i])->nnonz / (SCIP_Real) lpdata->ncols);
      }
      SCIPsortDownRealInt(zerohalfcutpriorities, sortedzerohalfcuts, nzerohalfcuts);
    
      /* check orthogonality */
      for( si = 0; si < nzerohalfcuts; ++si )
      {
         cutdatai = zerohalfcuts[sortedzerohalfcuts[si]];
         if( cutdatai->cut == NULL || !cutdatai->addedtolp )
            continue;      
         hasminorthogonality = TRUE;
         for( sj = 0; hasminorthogonality && sj < si; ++sj )
         {
            cutdataj = zerohalfcuts[sortedzerohalfcuts[sj]];
            if( cutdataj->cut != NULL && cutdataj->addedtolp )
               if( SCIPisLT(scip, SCIProwGetOrthogonality(cutdatai->cut , cutdataj->cut, ORTHOFUNC), 
                     MINORTHO) )
                  hasminorthogonality = FALSE;
         }

         /* add cut to LP */
         if( hasminorthogonality && cutdatai->addedtolp )
         {        
            SCIP_CALL(SCIPaddCut(scip, NULL, cutdatai->cut, sepadata->forcecutstolp));
            if( !cutdatai->islocal )
            {
               SCIP_CALL(SCIPaddPoolCut(scip, cutdatai->cut));
            }
            cutdatai->addedtolp = TRUE;
         }        
         else
         {
            cutdatai->addedtolp = FALSE;        
            nignoredcuts++;
         }
      } 
      nsepacuts -= nignoredcuts;

#ifdef ZEROHALF__PRINT_STATISTICS
      printZerohalfCutsStatistics(scip, sepadata, zerohalfcuts, nzerohalfcuts, sortedzerohalfcuts,
         zerohalfcutpriorities, NULL, nsepacuts);
#endif
#else /* new cutpool version: does not seem to be better */
      {
         int ncutpool;
         int ncutpoolold;

         ncutpool = 0;

         /* calculate score for all zerohalfcuts in cutpool */
         for( i = 0 ; i < nzerohalfcuts; ++i )
         {
            cutdatai = zerohalfcuts[i];

            if( cutdatai->cut == NULL || !cutdatai->addedtolp )
               continue;      

            sortedzerohalfcuts[ncutpool] = i;
            zerohalfcutminortho[ncutpool] = 1.0;
            zerohalfcutpriorities[ncutpool] = SCIPfloor(scip, violationbucketsize * cutdatai->efficacy)
               + (1.0 - (SCIP_Real) cutdatai->nnonz / (SCIP_Real) lpdata->ncols) 
               + SCIPfloor(scip, violationbucketsize * zerohalfcutminortho[ncutpool]);

            ncutpool++;
         }
         ncutpoolold = ncutpool;

         /* cut selection loop */
         while( ncutpool > 0 )
         {
            SCIP_Real bestscore;
            int bestpos;
            SCIP_Real priotmp;
            SCIP_Real minorthotmp;
            int sortidxtmp;

            assert(cutdatai->addedtolp);

            bestscore = SCIP_REAL_MIN;
            bestpos = -1;

            /* find cut with best (highest) score */
            for( i = 0; i < ncutpool; i++ )
            {
               /* check if cut is current best cut */
               if( zerohalfcutpriorities[i] > bestscore )
               {
                  bestscore = zerohalfcutpriorities[i];
                  bestpos = i;
               }
            }

            cutdatai = zerohalfcuts[sortedzerohalfcuts[bestpos]];

            /* add best cut to LP */
            SCIP_CALL(SCIPaddCut(scip, NULL, cutdatai->cut, sepadata->forcecutstolp));
            if( !cutdatai->islocal )
            {
               SCIP_CALL(SCIPaddPoolCut(scip, cutdatai->cut));
            }
            
            priotmp = zerohalfcutpriorities[bestpos];
            minorthotmp = zerohalfcutminortho[bestpos];
            sortidxtmp = sortedzerohalfcuts[bestpos];

            /* delete best cut from cutpool */
            zerohalfcutpriorities[bestpos] = zerohalfcutpriorities[ncutpool-1];
            zerohalfcutminortho[bestpos] = zerohalfcutminortho[ncutpool-1];
            sortedzerohalfcuts[bestpos] = sortedzerohalfcuts[ncutpool-1];
            ncutpool--;

            /* save data for statistic output */
            zerohalfcutpriorities[ncutpool] = priotmp;
            zerohalfcutminortho[ncutpool] = minorthotmp;
            sortedzerohalfcuts[ncutpool] = sortidxtmp;
            
            /* update orthogonalities of remaining cuts in cutpool */
            j = 0;
            while( j < ncutpool )
            {
               SCIP_Real thisortho;
               
               cutdataj = zerohalfcuts[sortedzerohalfcuts[j]];
               thisortho = SCIProwGetOrthogonality(cutdatai->cut , cutdataj->cut, ORTHOFUNC);

               if( thisortho < MINORTHO )
               {
                  priotmp = zerohalfcutpriorities[j];
                  minorthotmp = zerohalfcutminortho[j];
                  sortidxtmp = sortedzerohalfcuts[j];
 
                  /* delete cut from cutpool */
                  zerohalfcutpriorities[j] = zerohalfcutpriorities[ncutpool-1];
                  zerohalfcutminortho[j] = zerohalfcutminortho[ncutpool-1];
                  sortedzerohalfcuts[j] = sortedzerohalfcuts[ncutpool-1];
                  cutdataj->addedtolp = FALSE;        
                  ncutpool--;
                  nignoredcuts++;

                  /* save data for statistic output */
                  zerohalfcutpriorities[ncutpool] = priotmp;
                  zerohalfcutminortho[ncutpool] = minorthotmp;
                  sortedzerohalfcuts[ncutpool] = sortidxtmp;
               }
               else
               {
                  /* update score */
                  if( thisortho < zerohalfcutminortho[j] )
                  {
                     zerohalfcutminortho[j] = thisortho;
                     zerohalfcutpriorities[j] = SCIPfloor(scip, violationbucketsize * cutdataj->efficacy)
                        + (1.0 - (SCIP_Real) cutdataj->nnonz / (SCIP_Real) lpdata->ncols) 
                        + SCIPfloor(scip, violationbucketsize * zerohalfcutminortho[j]);
                  }
                  j++;
               }
            }
         }
         nsepacuts -= nignoredcuts;

#ifdef ZEROHALF__PRINT_STATISTICS
         printZerohalfCutsStatistics(scip, sepadata, zerohalfcuts, ncutpoolold, sortedzerohalfcuts,
            zerohalfcutpriorities, zerohalfcutminortho, nsepacuts);
#endif
      }      
#endif
    
      /* free temporary memory */
      SCIPfreeBufferArray(scip, &zerohalfcutminortho);
      SCIPfreeBufferArray(scip, &zerohalfcutpriorities);
      SCIPfreeBufferArray(scip, &sortedzerohalfcuts);
   } 
    
   sepadata->totalncutsfound += nzerohalfcuts;
   sepadata->totalnsepacuts += nsepacuts;

   /* free temporary memory */
   if( varsolvals != NULL )
   {
      SCIPfreeMemoryArray(scip, &varsolvals);
   }

   for( i = 0 ; i < nzerohalfcuts ; ++i)
   {
      SCIP_CALL( ZerohalfCutDataFree(scip, &(zerohalfcuts[i])) );  
   }
   SCIPfreeBufferArray(scip, &zerohalfcuts);
   SCIPfreeBufferArray(scip, &subproblempriorities);
   SCIPfreeBufferArray(scip, &sortedsubproblems);
   SCIP_CALL( ZerohalfLPDataFree(scip, &lpdata) );
  
   return SCIP_OKAY;
}


/* --------------------------------------------------------------------------------------------------------------------
 * separator specific interface methods 
 * -------------------------------------------------------------------------------------------------------------------- */

/** creates the zerohalf separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaZerohalf(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA*        sepadata;
   SCIP_SEPA* sepa;

   /* description of the preprocessing methods parameter */
   char preprocessingmethodsdescription[SCIP_MAXSTRLEN];
   /* description of the sepamethods parameter */
   char sepamethodsdescription[SCIP_MAXSTRLEN];
   /* description of the subscip parameter */
   char subscipobjectivedescription[SCIP_MAXSTRLEN];
   int ncharsprinted;

   ncharsprinted = SCIPmemccpy(preprocessingmethodsdescription,
      "preprocessing methods and ordering:\n"
      "   #                      'd' columns with small LP solution,\n"
      "   #                      'G' modified Gaussian elimination,\n"
      "   #                      'i' identical columns,\n"
      "   #                      'I' identical rows,\n"
      "   #                      'L' large slack rows,\n"
      "   #                      'M' large slack rows (minslack),\n"
      "   #                      's' column singletons,\n", '\0', SCIP_MAXSTRLEN - 1);

   assert(ncharsprinted > 0 && ncharsprinted < SCIP_MAXSTRLEN);

   ncharsprinted += SCIPmemccpy(&(preprocessingmethodsdescription[ncharsprinted - 1]),
      "   #                      'X' add trivial zerohalf cuts,\n"
      "   #                      'z' zero columns,\n"
      "   #                      'Z' zero rows,\n"
      "   #                      'C' fast {'z','s'},\n"
      "   #                      'R' fast {'Z','L','I'}\n"
      "   #\n"
      "   #                      '-' no preprocessing\n"
      "   #", '\0', (unsigned int) (SCIP_MAXSTRLEN - ncharsprinted - 1));

   assert(ncharsprinted > 0 && ncharsprinted < SCIP_MAXSTRLEN);

   ncharsprinted = SCIPmemccpy(sepamethodsdescription,
      "separating methods and ordering:\n"
      "   #                      '!' stop further processing if a cut was found,\n"
      "   #                      '2' exact polynomial time algorithm (only if matrix has max 2 odd entries per row),\n"
      "   #                      'e' enumeration heuristics (k=1: try all preprocessed rows),\n"
      "   #                      'E' enumeration heuristics (k=2: try all combinations of up to two preprocessed rows),\n"
      "   #                      'g' Extended Gaussian elimination heuristics,\n", '\0', SCIP_MAXSTRLEN - 1);

   assert(ncharsprinted > 0 && ncharsprinted < SCIP_MAXSTRLEN);

   ncharsprinted += SCIPmemccpy(&(sepamethodsdescription[ncharsprinted - 1]),
      "   #                      's' auxiliary IP heuristics (i.e. number of solved nodes is limited)\n"
      "   #                      'S' auxiliary IP exact      (i.e. unlimited number of nodes)\n"
      "   #\n"
      "   #                      '-' no processing\n"
      "   #", '\0', (unsigned int) (SCIP_MAXSTRLEN - ncharsprinted - 1));

   assert(ncharsprinted > 0 && ncharsprinted < SCIP_MAXSTRLEN);

   ncharsprinted = SCIPmemccpy(subscipobjectivedescription,
      "auxiliary IP objective:\n"
      "   #                      'v' maximize cut violation,\n"
      "   #                      'u' minimize number of aggregated rows in cut,\n"
      "   #                      'w' minimize number of aggregated rows in cut\n"
      "   #                          weighted by the number of rows in the aggregation,\n", '\0', SCIP_MAXSTRLEN - 1);

   assert(ncharsprinted > 0 && ncharsprinted < SCIP_MAXSTRLEN);

   ncharsprinted += SCIPmemccpy(&(subscipobjectivedescription[ncharsprinted - 1]),
      "   #                      'p' maximize cut violation and penalize a high number\n"
      "   #                          of aggregated rows in the cut weighted by the number\n"
      "   #                          of rows in the aggregation and the penalty factor p\n"
      "   #", '\0', (unsigned int) (SCIP_MAXSTRLEN - ncharsprinted - 1));

   assert(ncharsprinted > 0 && ncharsprinted < SCIP_MAXSTRLEN);

   /* create zerohalf separator data */
   SCIP_CALL(SCIPallocMemory(scip, &sepadata));

   /* initialize statistics */
   sepadata->totalncutsfound = 0;
   sepadata->totalnsepacuts = 0;
   sepadata->pptimers = NULL;
   sepadata->dtimer = NULL;
   sepadata->sepatimers = NULL;
   sepadata->nsepacutsalgo = NULL;
   sepadata->nzerohalfcutsalgo = NULL;
  
   sepadata->ppmethods = NULL;
   sepadata->sepamethods = NULL;
   sepadata->nppmethods = -1;
   sepadata->nsepamethods = -1;  
   sepadata->subscipsettings = NULL;

   sepadata->norigrows = 0;
   sepadata->origrows = NULL;
  
   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpZerohalf, NULL,
         sepadata) );

   assert(sepa != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyZerohalf) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeZerohalf) );
  
   /* add zerohalf separator parameters */
   SCIP_CALL(SCIPaddIntParam(scip,
         "separating/zerohalf/maxrounds",
         "maximal number of zerohalf separation rounds per node (-1: unlimited)",
         &(sepadata->maxrounds), FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip,
         "separating/zerohalf/maxroundsroot",
         "maximal number of zerohalf separation rounds in the root node (-1: unlimited)",
         &(sepadata->maxroundsroot), FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip,
         "separating/zerohalf/maxsepacuts",
         "maximal number of {0,1/2}-cuts separated per separation round",
         &(sepadata->maxsepacuts), FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip,
         "separating/zerohalf/maxsepacutsroot",
         "maximal number of {0,1/2}-cuts separated per separation round in the root node",
         &(sepadata->maxsepacutsroot), FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL));
   SCIP_CALL(SCIPaddBoolParam(scip,
         "separating/zerohalf/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &(sepadata->dynamiccuts), FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL));
  
  
   SCIP_CALL(SCIPaddIntParam(scip,
         "separating/zerohalf/maxcutsfound",
         "maximal number of {0,1/2}-cuts determined per separation round\n\
   #                      (this includes separated but inefficacious cuts)",
         &(sepadata->maxcuts), TRUE, DEFAULT_MAXCUTS, 0, 1e6, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip,
         "separating/zerohalf/maxcutsfoundroot",
         "maximal number of {0,1/2}-cuts determined per separation round in the root node\n\
   #                      (this includes separated but inefficacious cuts)",
         &(sepadata->maxcutsroot), TRUE, DEFAULT_MAXCUTSROOT, 0, 1e6, NULL, NULL));


   SCIP_CALL(SCIPaddIntParam(scip,
         "separating/zerohalf/maxdepth",
         "separating cuts only if depth <= maxdepth (-1: unlimited)",
         &(sepadata->maxdepth), TRUE, DEFAULT_MAXDEPTH, -1, INT_MAX, NULL, NULL));
   SCIP_CALL(SCIPaddLongintParam (scip,
         "separating/zerohalf/maxncalls", 
         "maximal number of calls (-1: unlimited)",
         &(sepadata->maxncalls), TRUE, DEFAULT_MAXNCALLS, -1LL, SCIP_LONGINT_MAX, NULL, NULL));
   SCIP_CALL(SCIPaddBoolParam(scip,
         "separating/zerohalf/relaxcontvars", 
         "should continuous variables be relaxed by adding variable bounds?", 
         &(sepadata->relaxcontvars), TRUE, DEFAULT_RELAXCONTVARS, NULL, NULL)); /**@todo support this */
   SCIP_CALL(SCIPaddBoolParam(scip,
         "separating/zerohalf/scalefraccoeffs",
         "should rows be scaled to make fractional coefficients integer?",  
         &(sepadata->scalefraccoeffs), TRUE, DEFAULT_SCALEFRACCOEFFS, NULL, NULL));
   SCIP_CALL(SCIPaddBoolParam(scip,
         "separating/zerohalf/ignoreprevzhcuts",
         "should zerohalf cuts found in previous callbacks ignored?",
         &(sepadata->ignoreprevzhcuts), TRUE, DEFAULT_IGNOREPREVIOUSZHCUTS, NULL, NULL));
   SCIP_CALL(SCIPaddBoolParam(scip,
         "separating/zerohalf/onlyorigrows",
         "should only original LP rows be considered (i.e. ignore previously added LP rows)?",
         &(sepadata->onlyorigrows), TRUE, DEFAULT_ONLYORIGROWS, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
         "separating/zerohalf/usezhcutpool",
         "should zerohalf cuts be filtered using a cutpool?",
         &(sepadata->usezhcutpool), TRUE, DEFAULT_USEZHCUTPOOL, NULL, NULL));
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/zerohalf/maxtestdelta",
         "maximal number of different deltas to try for cmir (-1: unlimited, 0: delta=1)",
         &sepadata->maxtestdelta, TRUE, DEFAULT_MAXTESTDELTA, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/zerohalf/trynegscaling",
         "should negative values also be tested in scaling for cmir?",
         &sepadata->trynegscaling, TRUE, DEFAULT_TRYNEGSCALING, NULL, NULL) );
 
  
   SCIP_CALL(SCIPaddBoolParam(scip,
         "separating/zerohalf/preprocessing/decomposeproblem",
         "should problem be decomposed into subproblems (if possible) before applying preprocessing?",
         &(sepadata->decomposeproblem), FALSE, DEFAULT_DECOMPOSEPROBLEM, NULL, NULL)); 
   SCIP_CALL(SCIPaddRealParam(scip,
         "separating/zerohalf/preprocessing/delta",
         "value of delta parameter used in preprocessing method 'd'",
         &(sepadata->ppdelta), FALSE, DEFAULT_PPDELTA, 0.0, 1.0, NULL, NULL));
   SCIP_CALL(SCIPaddStringParam(scip,
         "separating/zerohalf/preprocessing/ppmethods",
         preprocessingmethodsdescription,
         &(sepadata->ppmethods), FALSE, DEFAULT_PPMETHODS, NULL, NULL));


  
   SCIP_CALL(SCIPaddBoolParam(scip,
         "separating/zerohalf/separating/forcecutstolp",
         "should the cuts be forced to enter the LP?",
         &(sepadata->forcecutstolp), FALSE, DEFAULT_FORCECUTSTOLP, NULL, NULL));
   SCIP_CALL(SCIPaddBoolParam(scip,
         "separating/zerohalf/separating/forcecutstosepastore",
         "should the cuts be forced to enter SCIP's sepastore?",
         &(sepadata->forcecutstosepastore), FALSE, DEFAULT_FORCECUTSTOSEPASTORE, NULL, NULL));
   SCIP_CALL(SCIPaddRealParam(scip,
         "separating/zerohalf/separating/minviolation",
         "minimal violation of a {0,1/2}-cut to be separated",
         &(sepadata->minviolation), FALSE, DEFAULT_MINVIOLATION, 0.001, 0.5, NULL, NULL));
   SCIP_CALL(SCIPaddStringParam(scip,
         "separating/zerohalf/separating/sepamethods",
         sepamethodsdescription,
         &(sepadata->sepamethods), FALSE, DEFAULT_SEPAMETHODS, NULL, NULL));


  
   SCIP_CALL(SCIPaddStringParam(scip,
         "separating/zerohalf/separating/auxip/settingsfile",
         "optional settings file of the auxiliary IP (-: none)",
         &(sepadata->subscipsettings), FALSE, DEFAULT_SUBSCIPSETTINGS, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip,
         "separating/zerohalf/separating/auxip/sollimit",
         "limits/solutions setting of the auxiliary IP",
         &(sepadata->subscipsollimit), FALSE, DEFAULT_SUBSCIPSOLLIMIT, -1, INT_MAX, NULL, NULL));
   SCIP_CALL(SCIPaddRealParam(scip,
         "separating/zerohalf/separating/auxip/penaltyfactor",
         "penalty factor used with objective function 'p' of auxiliary IP",
         &(sepadata->subscipobjpen), FALSE, DEFAULT_SUBSCIPOBJPEN, 0.0, 1.0, NULL, NULL));
   SCIP_CALL(SCIPaddBoolParam(scip,
         "separating/zerohalf/separating/auxip/useallsols", 
         "should all (proper) solutions of the auxiliary IP be used to generate cuts instead of using only the best?",
         &(sepadata->subscipuseallsols), FALSE, DEFAULT_SUBSCIPUSEALLSOLS, NULL, NULL));
   SCIP_CALL(SCIPaddCharParam(scip,
         "separating/zerohalf/separating/auxip/objective", 
	 subscipobjectivedescription,
         &(sepadata->subscipobjective), FALSE, DEFAULT_SUBSCIPOBJECTIVE, "uvwp", NULL, NULL));

   return SCIP_OKAY;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_bitvar.c
 * @brief  constraint handler for bitvar constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_bitvar.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "bitvar"
#define CONSHDLR_DESC          "arbitrarily long integer variables represented as bit strings"
#define CONSHDLR_SEPAPRIORITY  +2000000
#define CONSHDLR_ENFOPRIORITY  - 500000
#define CONSHDLR_CHECKPRIORITY - 500000
#define CONSHDLR_SEPAFREQ             1
#define CONSHDLR_PROPFREQ             1
#define CONSHDLR_NEEDSCONS         TRUE

#define EVENTHDLR_NAME         "bitvar"
#define EVENTHDLR_DESC         "bound change event handler for bitvar constraints"

#define WORDSIZE                      8 /**< number of bits in one word of the bitvar */
#define WORDPOWER       (1 << WORDSIZE) /**< number of different values of one word (2^WORDSIZE) */



/*
 * Data structures
 */

/** constraint data for bitvar constraints */
struct ConsData
{
   VAR**            bits;               /**< binaries representing bits of the bitvar, least significant first */
   VAR**            words;              /**< integers representing words of the bitvar, least significant first */
   ROW**            rows;               /**< LP rows storing equalities for each word */
   int              nbits;              /**< number of bits */
   int              nwords;             /**< number of words: nwords = ceil(nbits/WORDSIZE) */
   unsigned int     constant:1;         /**< was bitvar created as a constant? */
   unsigned int     propagated:1;       /**< is constraint already preprocessed/propagated? */
};

/** constraint handler data */
struct ConsHdlrData
{
   EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
};




/*
 * Local methods
 */

/** returns the number of bits of the given word */
static
int wordSize(
   CONSDATA*        consdata,           /**< bitvar constraint data */
   int              word                /**< word number */
   )
{
   assert(consdata != NULL);
   assert(0 <= word && word < consdata->nwords);

   if( word < consdata->nwords-1 )
      return WORDSIZE;
   else
      return consdata->nbits - (consdata->nwords-1) * WORDSIZE;
}

/** returns the number of different values the given word can store (2^#bits) */
static
int wordPower(
   CONSDATA*        consdata,           /**< bitvar constraint data */
   int              word                /**< word number */
   )
{
   assert(consdata != NULL);
   assert(0 <= word && word < consdata->nwords);

   if( word < consdata->nwords-1 )
      return WORDPOWER;
   else
      return (1 << (consdata->nbits - (consdata->nwords-1) * WORDSIZE));
}

/** creates constaint handler data for bitvar constraint handler */
static
RETCODE conshdlrdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   CHECK_OKAY( SCIPallocMemory(scip, conshdlrdata) );

   /* get event handler for processing bound change events */
   (*conshdlrdata)->eventhdlr = SCIPfindEventHdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      errorMessage("event handler for bitvar constraints not found");
      return SCIP_PLUGINNOTFOUND;
   }

   return SCIP_OKAY;
}

/** frees constaint handler data for bitvar constraint handler */
static
RETCODE conshdlrdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSHDLRDATA**   conshdlrdata        /**< pointer the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   SCIPfreeMemory(scip, conshdlrdata);

   return SCIP_OKAY;
}

/** creates a bitvar constraint data object along with the corresponding variables */
static
RETCODE consdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to store the bitvar constraint data */
   int              nbits               /**< number of bits in the bitvar */
   )
{
   assert(consdata != NULL);
   assert(nbits >= 1);

   /* allocate memory */
   CHECK_OKAY( SCIPallocBlockMemory(scip, consdata) );
   (*consdata)->nbits = nbits;
   (*consdata)->nwords = (nbits+WORDSIZE-1)/WORDSIZE;
   CHECK_OKAY( SCIPallocBlockMemoryArray(scip, &(*consdata)->bits, (*consdata)->nbits) );
   CHECK_OKAY( SCIPallocBlockMemoryArray(scip, &(*consdata)->words, (*consdata)->nwords) );
   (*consdata)->rows = NULL;
   (*consdata)->constant = FALSE;
   (*consdata)->propagated = FALSE;

   /* clear the bits and words arrays */
   clearMemoryArray((*consdata)->bits, (*consdata)->nbits);
   clearMemoryArray((*consdata)->words, (*consdata)->nwords);

   return SCIP_OKAY;
}

/** creates variables for the bitvar and adds them to the problem */
static
RETCODE consdataCreateVars(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< bitvar constraint data */
   EVENTHDLR*       eventhdlr,          /**< event handler for bound change events */
   const char*      name,               /**< prefix for variable names */
   Real             obj                 /**< objective value of bitvar variable */
   )
{
   char varname[MAXSTRLEN];
   Real bitobj;
   int i;

   assert(consdata != NULL);
   assert(name != NULL);

   /* create binary variables for bits */
   bitobj = obj;
   for( i = 0; i < consdata->nbits; ++i )
   {
      sprintf(varname, "%s_b%d", name, i);
      CHECK_OKAY( SCIPcreateVar(scip, &consdata->bits[i], varname, 0.0, 1.0, bitobj, SCIP_VARTYPE_BINARY, TRUE) );
      CHECK_OKAY( SCIPaddVar(scip, consdata->bits[i]) );
      bitobj *= 2.0;
      
      debugMessage("created bit variable <%s> with obj=%g for bitvar constraint\n",
         SCIPvarGetName(consdata->bits[i]), bitobj);

      /* if we are in the transformed problem, catch bound tighten events on variable */
      if( SCIPvarIsTransformed(consdata->bits[i]) )
      {
         CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->bits[i], SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr,
                        (EVENTDATA*)consdata) );
      }
   }

   /* create integer variables for words */
   for( i = 0; i < consdata->nwords; ++i )
   {
      sprintf(varname, "%s_w%d", name, i);
      CHECK_OKAY( SCIPcreateVar(scip, &consdata->words[i], varname, 0.0, wordPower(consdata, i)-1.0, 0.0,
                     SCIP_VARTYPE_INTEGER, TRUE) );
      CHECK_OKAY( SCIPaddVar(scip, consdata->words[i]) );
      
      debugMessage("created word variable <%s> for bitvar constraint\n", SCIPvarGetName(consdata->words[i]));

      /* if we are in the transformed problem, catch bound tighten events on variable */
      if( SCIPvarIsTransformed(consdata->words[i]) )
      {
         CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->words[i], SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr,
                        (EVENTDATA*)consdata) );
      }
   }

   /* issue warning message, if objective value for most significant bit grew too large */
   if( ABS(bitobj) >= SCIPinfinity(scip) )
   {
      char msg[MAXSTRLEN];

      sprintf(msg, "objective value %g of %d-bit variable grew up to %g in last bit (exceeds infinity)\n", 
         obj, consdata->nbits, bitobj/2.0);
      errorMessage(msg);
      return SCIP_INVALIDDATA;
   }
   if( ABS(bitobj) > SCIPinfinity(scip)/10000.0 )
   {
      char msg[MAXSTRLEN];

      sprintf(msg, "Warning! objective value %g of %d-bit variable grew up to %g in last bit\n", 
         obj, consdata->nbits, bitobj/2.0);
      SCIPmessage(scip, SCIP_VERBLEVEL_MINIMAL, msg);
   }
   
   return SCIP_OKAY;
}

/** create variables in target constraint data by transforming variables of the source constraint data */
static
RETCODE consdataTransformVars(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        sourcedata,         /**< bitvar constraint data of the source constraint */
   CONSDATA*        targetdata,         /**< bitvar constraint data of the target constraint */
   EVENTHDLR*       eventhdlr           /**< event handler for bound change events */
   )
{
   int i;

   assert(sourcedata != NULL);
   assert(sourcedata->bits != NULL);
   assert(sourcedata->words != NULL);
   assert(targetdata != NULL);
   assert(targetdata->bits != NULL);
   assert(targetdata->words != NULL);
   assert(sourcedata->nbits == targetdata->nbits);
   assert(sourcedata->nwords == targetdata->nwords);

   /* get transformed variables */
   CHECK_OKAY( SCIPgetTransformedVars(scip, sourcedata->nbits, sourcedata->bits, targetdata->bits) );
   CHECK_OKAY( SCIPgetTransformedVars(scip, sourcedata->nwords, sourcedata->words, targetdata->words) );

   /* capture bit variables and catch bound tighten events */
   for( i = 0; i < sourcedata->nbits; ++i )
   {
      assert(SCIPvarIsTransformed(targetdata->bits[i]));
      CHECK_OKAY( SCIPcaptureVar(scip, targetdata->bits[i]) );
      CHECK_OKAY( SCIPcatchVarEvent(scip, targetdata->bits[i], SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr,
                     (EVENTDATA*)targetdata) );
   }

   /* capture word variables and catch bound tighten events */
   for( i = 0; i < sourcedata->nwords; ++i )
   {
      assert(SCIPvarIsTransformed(targetdata->words[i]));
      CHECK_OKAY( SCIPcaptureVar(scip, targetdata->words[i]) );
      CHECK_OKAY( SCIPcatchVarEvent(scip, targetdata->words[i], SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr,
                     (EVENTDATA*)targetdata) );
   }

   return SCIP_OKAY;
}

/** frees a bitvar constraint data object and releases corresponding variables */
static
RETCODE consdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to the bitvar constraint data */
   EVENTHDLR*       eventhdlr           /**< event handler for bound change events */
   )
{
   int i;

   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release binary variables for bits */
   for( i = 0; i < (*consdata)->nbits; ++i )
   {
      /* if we are in the transformed problem, drop bound tighten events on variable */
      if( SCIPvarIsTransformed((*consdata)->bits[i]) )
      {
         CHECK_OKAY( SCIPdropVarEvent(scip, (*consdata)->bits[i], eventhdlr, (EVENTDATA*)(*consdata)) );
      }

      /* release variable */
      CHECK_OKAY( SCIPreleaseVar(scip, &(*consdata)->bits[i]) );
   }

   /* release integer variables for words */
   for( i = 0; i < (*consdata)->nwords; ++i )
   {
      /* if we are in the transformed problem, drop bound tighten events on variable */
      if( SCIPvarIsTransformed((*consdata)->words[i]) )
      {
         CHECK_OKAY( SCIPdropVarEvent(scip, (*consdata)->words[i], eventhdlr, (EVENTDATA*)(*consdata)) );
      }

      /* release variable */
      CHECK_OKAY( SCIPreleaseVar(scip, &(*consdata)->words[i]) );
   }

   /* release rows */
   if( (*consdata)->rows != NULL )
   {
      for( i = 0; i < (*consdata)->nwords; ++i )
      {
         if( (*consdata)->rows[i] != NULL )
         {
            CHECK_OKAY( SCIPreleaseRow(scip, &(*consdata)->rows[i]) );
         }
      }
   }

   /* free memory */
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->bits, (*consdata)->nbits);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->words, (*consdata)->nwords);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->rows, (*consdata)->nwords);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** checks given word of bitvar constraint for feasibility of given solution or actual solution */
static
RETCODE checkWord(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitvar constraint */
   int              word,               /**< word number to check */
   SOL*             sol,                /**< solution to be checked, or NULL for actual solution */
   Bool             checklprows,        /**< has bitvar constraint to be checked, if it is already in current LP? */
   int*             nviolatedbits       /**< pointer to store the number of violated bits */
   )
{
   CONSDATA* consdata;
   Real wordsol;
   Real bitsol;
   Bool wordbitisset;
   Bool bitsolisone;
   int wordsolint;
   int wordsize;
   int bitmask;
   int b;

   assert(nviolatedbits != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= word && word < consdata->nwords);

   debugMessage("checking bitvar constraint <%s> at word %d\n", SCIPconsGetName(cons), word);

   *nviolatedbits = 0;

   if( !checklprows && consdata->rows != NULL && consdata->rows[word] != NULL && SCIProwIsInLP(consdata->rows[word]) )
      return SCIP_OKAY;

   /* get the value of the word and convert it into an integer */
   wordsol = SCIPgetSolVal(scip, sol, consdata->words[word]);
   assert(SCIPisIntegral(scip, wordsol));
   wordsolint = (int)(wordsol+0.5);
   assert(SCIPisFeasEQ(scip, wordsol, (Real)wordsolint));

   /* compare each bit in the word's solution with the value of the corresponding binary variable */
   wordsize = wordSize(consdata, word);
   bitmask = 1;
   for( b = 0; b < wordsize; ++b )
   {
      assert(0 <= bitmask && bitmask <= WORDPOWER/2);
      bitsol = SCIPgetSolVal(scip, sol, consdata->bits[word*WORDSIZE+b]);
      assert(SCIPisIntegral(scip, bitsol));
      assert(SCIPisFeasEQ(scip, bitsol, 0.0) || SCIPisFeasEQ(scip, bitsol, 1.0));
      bitsolisone = (bitsol > 0.5);
      wordbitisset = ((wordsolint & bitmask) > 0);
      if( bitsolisone != wordbitisset )
         (*nviolatedbits)++;
      bitmask <<= 1;
   }

   /* update constraint's age */
   if( *nviolatedbits == 0 )
   {
      CHECK_OKAY( SCIPincConsAge(scip, cons) );
   }
   else
   {
      CHECK_OKAY( SCIPresetConsAge(scip, cons) );
   }

   return SCIP_OKAY;
}

/** checks all words of bitvar constraint for feasibility of given solution or actual solution */
static
RETCODE checkCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitvar constraint */
   SOL*             sol,                /**< solution to be checked, or NULL for actual solution */
   Bool             checklprows,        /**< has bitvar constraint to be checked, if it is already in current LP? */
   Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   CONSDATA* consdata;
   int nviolatedbits;
   int w;

   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violated = FALSE;

   for( w = 0; w < consdata->nwords && !(*violated); ++w )
   {
      CHECK_OKAY( checkWord(scip, cons, w, sol, checklprows, &nviolatedbits) );
      (*violated) |= (nviolatedbits > 0);
   }

   return SCIP_OKAY;
}

/** creates an LP row for a single word in a bitvar constraint */
static
RETCODE createRow(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitvar constraint */
   int              word                /**< word number to create the row for */
   )
{
   CONSDATA* consdata;
   char rowname[MAXSTRLEN];
   Real coef;
   int bitstart;
   int bitend;
   int b;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->bits != NULL);
   assert(consdata->words != NULL);
   assert(0 <= word && word < consdata->nwords);

   /* create the rows array, if not yet existing */
   if( consdata->rows == NULL )
   {
      CHECK_OKAY( SCIPallocBlockMemoryArray(scip, &consdata->rows, consdata->nwords) );
      clearMemoryArray(consdata->rows, consdata->nwords);
   }
   assert(consdata->rows != NULL);
   assert(consdata->rows[word] == NULL);

   /* create equality  - word + 2^0*bit[0] + 2^1*bit[1] + ... + 2^(WORDSIZE-1)*bit[WORDSIZE-1] == 0 */
   sprintf(rowname, "c_%s", SCIPvarGetName(consdata->words[word]));
   CHECK_OKAY( SCIPcreateEmptyRow(scip, &consdata->rows[word], rowname, 0.0, 0.0,
                  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   
   CHECK_OKAY( SCIPaddVarToRow(scip, consdata->rows[word], consdata->words[word], -1.0) );
   bitstart = word*WORDSIZE;
   bitend = bitstart + wordSize(consdata, word);
   assert(bitstart < bitend);
   coef = 1.0;
   for( b = bitstart; b < bitend; ++b )
   {
      CHECK_OKAY( SCIPaddVarToRow(scip, consdata->rows[word], consdata->bits[b], coef) );
      coef *= 2.0;
   }
   assert(SCIPisEQ(scip, coef, (Real)(wordPower(consdata, word))));

   return SCIP_OKAY;
}

/** adds bitvar constraint as cuts to the LP */
static
RETCODE addCut(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitvar constraint */
   int              word,               /**< word number to add the cut for */
   Real             violation           /**< absolute violation of the constraint */
   )
{
   CONSDATA* consdata;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= word && word < consdata->nwords);

   /* create the selected row, if not yet existing */
   if( consdata->rows == NULL || consdata->rows[word] == NULL )
   {
      /* convert consdata object into LP row */
      CHECK_OKAY( createRow(scip, cons, word) );
   }
   assert(consdata->rows != NULL);
   assert(consdata->rows[word] != NULL);
   assert(!SCIProwIsInLP(consdata->rows[word]));
   
   /* insert LP row as cut */
   CHECK_OKAY( SCIPaddCut(scip, consdata->rows[word], 
                  violation/SCIProwGetNorm(consdata->rows[word])/(SCIProwGetNNonz(consdata->rows[word])+1)) );

   return SCIP_OKAY;
}

/** separates bitvar constraint: adds each word of bitvar constraint as cut, if violated by current LP solution */
static
RETCODE separateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitvar constraint */
   RESULT*          result              /**< pointer to store result of separation */
   )
{
   CONSDATA* consdata;
   int w;

   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( w = 0; w < consdata->nwords; ++w )
   {
      int nviolatedbits;

      CHECK_OKAY( checkWord(scip, cons, w, NULL, FALSE, &nviolatedbits) );
      
      if( nviolatedbits > 0 )
      {
         /* insert LP row as cut */
         CHECK_OKAY( addCut(scip, cons, w, (Real)nviolatedbits) );
         *result = SCIP_SEPARATED;
      }
   }

   return SCIP_OKAY;
}

/** propagates domains of variables of a single word in a bitvar constraint */
static
RETCODE propagateWord(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitvar constraint */
   int              word,               /**< word number to add the cut for */
   int*             nfixedvars,         /**< pointer to add up the number of fixed variables, or NULL */
   int*             nchgbds,            /**< pointer to add up the number of fixed bounds, or NULL */
   Bool*            infeasible,         /**< pointer to store whether the constraint is infeasible in current bounds */
   Bool*            wordfixed           /**< pointer to store whether the whole word is fixed */
   )
{
   CONSDATA* consdata;
   VAR* wordvar;
   Real lb;
   Real ub;
   Bool lbbit;
   Bool ubbit;
   int oldwordlb;
   int oldwordub;
   int wordlb;
   int wordub;
   int bitslb;
   int bitsub;
   int wordsize;
   int bitstart;
   int bitend;
   int bitval;
   int b;

   assert(infeasible != NULL);
   assert(wordfixed != NULL);

   *infeasible = FALSE;
   *wordfixed = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->bits != NULL);
   assert(consdata->words != NULL);

   /* get bit positions */
   wordsize = wordSize(consdata, word);
   bitstart = word*WORDSIZE;
   bitend = bitstart + wordsize;
   assert(wordsize >= 1);
   assert(bitstart < bitend);

   /* get the word along with its bounds */
   wordvar = consdata->words[word];
   lb = SCIPvarGetLbLocal(wordvar);
   ub = SCIPvarGetUbLocal(wordvar);
   assert(SCIPisIntegral(scip, lb));
   assert(SCIPisIntegral(scip, ub));
   wordlb = (int)(lb + 0.5);
   wordub = (int)(ub + 0.5);
   oldwordlb = wordlb;
   oldwordub = wordub;

   /* scan the bits for fixed bits, calculating the lowest and largest possible value with these fixings */
   bitslb = 0;
   bitsub = 0;
   for( b = bitend-1, bitval = (1 << (wordsize-1)); b >= bitstart; --b, bitval >>= 1 )
   {
      lb = SCIPvarGetLbLocal(consdata->bits[b]);
      ub = SCIPvarGetUbLocal(consdata->bits[b]);
      if( lb > 0.5 )
      {
         /* bit is fixed to 1: increase the lower bound appropriately */
         bitslb |= bitval;
      }
      if( ub > 0.5 )
      {
         /* bit is not fixed to 0: increase the upper bound appropriately */
         bitsub |= bitval;
      }
   }
   assert(bitval == 0);

   /* update the word's bounds */
   wordlb = MAX(wordlb, bitslb);
   wordub = MIN(wordub, bitsub);

   /* scan the bits from most to least significant bit (from right to left), and check for possible bit values */
   /*debugMessage("bitvar <%s> word %d: prop loop start: word=[%d,%d], bits=[%d,%d] -- [0x%x,0x%x] / [0x%x,0x%x]\n",
     SCIPconsGetName(cons), word, wordlb, wordub, bitslb, bitsub, wordlb, wordub, bitslb, bitsub);*/
   for( b = bitend-1, bitval = (1 << (wordsize-1));
        b >= bitstart && wordlb < wordub && (bitslb < wordlb || bitsub > wordub)
           && (bitsub - bitval < wordlb || bitslb + bitval > wordub);
        --b, bitval >>= 1 )
   {
      /*debugMessage("bitvar <%s> word %d: prop bit %d: word=[%d,%d], bits=[%d,%d] -- [0x%x,0x%x] / [0x%x,0x%x]\n",
        SCIPconsGetName(cons), word, b, wordlb, wordub, bitslb, bitsub, wordlb, wordub, bitslb, bitsub);*/

      /* check, if the largest possible value with the current bit set to 0 violates the word's lower bound */
      if( (bitsub & ~bitval) < wordlb )
      {
         /* the bit can't be set to 0: fix the bit to 1 in the bit's lower bound and update the word's lower bound */
         bitslb |= bitval;
         wordlb = MAX(wordlb, bitslb);
      }
      /* check, if the smallest possible value with the current bit set to 1 violates the word's upper bound */
      if( (bitslb | bitval) > wordub )
      {
         /* the bit can't be set to 1: fix the bit to 0 in the bit's upper bound and update the word's upper bound */
         bitsub &= ~bitval;
         wordub = MIN(wordub, bitsub);
      }
      assert(bitslb <= wordlb && bitsub >= wordub);
   }
   assert(wordlb >= oldwordlb);
   assert(wordub <= oldwordub);
   /*debugMessage("bitvar <%s> word %d: prop loop end: word=[%d,%d], bits=[%d,%d] -- [0x%x,0x%x] / [0x%x,0x%x]\n",
     SCIPconsGetName(cons), word, wordlb, wordub, bitslb, bitsub, wordlb, wordub, bitslb, bitsub);*/

   /* change the word's bounds */
   if( wordlb > wordub )
   {
      debugMessage("bitvar <%s> word %d: infeasibility at bit %d detected: word=[%d,%d] bits=[%d,%d]\n",
         SCIPconsGetName(cons), word, b, wordlb, wordub, bitslb, bitsub);
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   else if( wordlb == wordub )
   {
      if( wordlb > oldwordlb || wordub < oldwordub )
      {
         debugMessage("bitvar <%s> word %d: fixing word <%s> to %d\n",
            SCIPconsGetName(cons), word, SCIPvarGetName(wordvar), wordlb);
         CHECK_OKAY( SCIPfixVar(scip, wordvar, (Real)wordlb, infeasible) );
         if( *infeasible )
            return SCIP_OKAY;
         if( nfixedvars != NULL )
            (*nfixedvars)++;
      }
      *wordfixed = TRUE;
   }
   else
   {
      if( wordlb > oldwordlb )
      {
         debugMessage("bitvar <%s> word %d: tighten word <%s> lower bound: [%d,%d] -> [%d,%d]\n",
            SCIPconsGetName(cons), word, SCIPvarGetName(wordvar), oldwordlb, oldwordub, wordlb, oldwordub);
         CHECK_OKAY( SCIPchgVarLb(scip, wordvar, (Real)wordlb) );
         if( nchgbds != NULL )
            (*nchgbds)++;
      }
      if( wordub < oldwordub )
      {
         debugMessage("bitvar <%s> word %d: tighten word <%s> upper bound: [%d,%d] -> [%d,%d]\n",
            SCIPconsGetName(cons), word, SCIPvarGetName(wordvar), wordlb, oldwordub, wordlb, wordub);
         CHECK_OKAY( SCIPchgVarUb(scip, wordvar, (Real)wordub) );
         if( nchgbds != NULL )
            (*nchgbds)++;
      }
   }

   /* fix most significant bit variables as long as both bounds are equal on these bits */
   for( b = bitend-1, bitval = (1 << (wordsize-1)); b >= bitstart; --b, bitval >>= 1 )
   {
      lbbit = ((wordlb & bitval) > 0);
      ubbit = ((wordub & bitval) > 0);
      assert(!lbbit || ubbit); /* otherwise, wordlb > wordub */
      if( lbbit != ubbit )
         break; /* stop fixing at first bit discrepancy */
      if( !lbbit )
      {
         /* fix the bit variable to 0 */
         assert(SCIPvarGetLbLocal(consdata->bits[b]) < 0.5); /* otherwise, wordlb > wordub */
         if( SCIPvarGetUbLocal(consdata->bits[b]) > 0.5 )
         {
            debugMessage("bitvar <%s> word %d: fixing bit %d <%s> to 0 (word <%s>: [0x%x,0x%x])\n",
               SCIPconsGetName(cons), word, b, SCIPvarGetName(consdata->bits[b]),
               SCIPvarGetName(wordvar), wordlb, wordub);
            CHECK_OKAY( SCIPfixVar(scip, consdata->bits[b], 0.0, infeasible) );
            if( *infeasible )
               return SCIP_OKAY;
            if( nfixedvars != NULL )
               (*nfixedvars)++;
         }
      }
      else
      {
         /* fix the bit variable to 1 */
         assert(SCIPvarGetUbLocal(consdata->bits[b]) > 0.5); /* otherwise, wordlb > wordub */
         if( SCIPvarGetLbLocal(consdata->bits[b]) < 0.5 )
         {
            debugMessage("bitvar <%s> word %d: fixing bit %d <%s> to 1 (word <%s>: [0x%x,0x%x])\n",
               SCIPconsGetName(cons), word, b, SCIPvarGetName(consdata->bits[b]),
               SCIPvarGetName(wordvar), wordlb, wordub);
            CHECK_OKAY( SCIPfixVar(scip, consdata->bits[b], 1.0, infeasible) );
            if( *infeasible )
               return SCIP_OKAY;
            if( nfixedvars != NULL )
               (*nfixedvars)++;
         }
      }
   }
   assert((b < bitstart) == *wordfixed);

   return SCIP_OKAY;
}

/** propagates domains of variables of a bitvar constraint */
static
RETCODE propagateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitvar constraint */
   int*             nfixedvars,         /**< pointer to add up the number of fixed variables, or NULL */
   int*             nchgbds,            /**< pointer to add up the number of fixed bounds, or NULL */
   int*             ndelconss,          /**< pointer to add up the number of deleted constraints, or NULL */
   Bool*            infeasible          /**< pointer to store whether the constraint is infeasible in current bounds */
   )
{
   CONSDATA* consdata;
   Bool wordfixed;
   Bool allfixed;
   int w;

   assert(infeasible != NULL);

   *infeasible = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check, if the constraint is already propagated */
   if( consdata->propagated )
      return SCIP_OKAY;

   /* mark the constraint propagated */
   consdata->propagated = TRUE;

   /* propagate all words in the bitvar */
   allfixed = TRUE;
   for( w = 0; w < consdata->nwords && !(*infeasible); ++w )
   {
      CHECK_OKAY( propagateWord(scip, cons, w, nfixedvars, nchgbds, infeasible, &wordfixed) );
      allfixed &= wordfixed;
   }

   /* if all words in the bitvar are fixed, we don't need the constraint anymore */
   if( allfixed )
   {
      CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
      if( ndelconss != NULL )
         (*ndelconss)++;
   }

   return SCIP_OKAY;
}




/*
 * parsing strings to bitvar constants
 */

/** parses a string 'X...X' (X in {0,1}) with binary base into an array of bits, storing the
 *  bits with least significant bit first
 */
static
RETCODE parseBitString(
   const char*      s,                  /**< string to parse into bit array */
   int              nbits,              /**< number of bits to parse */
   Bool*            bitarray            /**< array to store the parsed bits in (must have length nbits) */
   )
{
   char msg[MAXSTRLEN];
   int b;
   int i;

   assert(s != NULL);
   assert(nbits >= 1);
   assert(bitarray != NULL);

   /* initialize the bitarray with FALSE's */
   clearMemoryArray(bitarray, nbits);

   /* scan the string from back to front */
   for( b = 0, i = strlen(s)-1; i >= 0; --i, ++b )
   {
      switch( s[i] )
      {
      case '0': /* nothing to do here: even allow additional 0's at the front */
         break;
      case '1':
         if( b >= nbits )
         {
            sprintf(msg, "bit constant too large for given bit size %d (digit %d sets bit %d)", nbits, i, b);
            errorMessage(msg);
            return SCIP_PARSEERROR;
         }
         bitarray[b] = TRUE;
         break;
      default:
         sprintf(msg, "invalid character <%c> at digit %d in binary string constant", s[i], i);
         errorMessage(msg);
         return SCIP_PARSEERROR;
      }
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
DECL_CONSFREE(consFreeBitvar)
{
   CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdataFree(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called when problem solving starts) */
#define consInitBitvar NULL


/** deinitialization method of constraint handler (called when problem solving exits) */
#define consExitBitvar NULL


/** solving start notification method of constraint handler (called when presolving was finished) */
#define consSolstartBitvar NULL


/** frees specific constraint data */
static
DECL_CONSDELETE(consDeleteBitvar)
{
   CONSHDLRDATA* conshdlrdata;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* free constraint data */
   CHECK_OKAY( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
DECL_CONSTRANS(consTransBitvar)
{
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* sourcedata;
   CONSDATA* targetdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* get source constraint data */
   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->rows == NULL);  /* in original problem, there cannot be LP rows */

   /* create constraint data for target constraint */
   CHECK_OKAY( consdataCreate(scip, &targetdata, sourcedata->nbits) );
   CHECK_OKAY( consdataTransformVars(scip, sourcedata, targetdata, conshdlrdata->eventhdlr) );

   /* create target constraint */
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                  SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
                  SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
                  SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), SCIPconsIsRemoveable(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
DECL_CONSINITLP(consInitlpBitvar)
{
   int c;
   int w;

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* add LP relaxations for all initial constraints */
   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsIsInitial(conss[c]) )
      {
         CONSDATA* consdata;

         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);

         /** add rows for all words in the bitvar constraint */
         debugMessage("adding initial bitvar constraint <%s> to LP\n", SCIPconsGetName(conss[c]));
         for( w = 0; w < consdata->nwords; ++w )
         {
            CHECK_OKAY( addCut(scip, conss[c], w, 0.0) );
         }
      }
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler */
static
DECL_CONSSEPA(consSepaBitvar)
{
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* step 1: check all useful bitvar constraints for feasibility */
   for( c = 0; c < nusefulconss; ++c )
   {
      debugMessage("separating bitvar constraint <%s>\n", SCIPconsGetName(conss[c]));
      CHECK_OKAY( separateCons(scip, conss[c], result) );
   }

   /* step 2: if no cuts were found and we are in the root node, check remaining bitvar constraints for feasibility */
   if( SCIPgetActDepth(scip) == 0 )
   {
      for( c = nusefulconss; c < nconss && *result == SCIP_DIDNOTFIND; ++c )
      {
         debugMessage("separating bitvar constraint <%s>\n", SCIPconsGetName(conss[c]));
         CHECK_OKAY( separateCons(scip, conss[c], result) );
      }
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpBitvar)
{
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* check for violated constraints
    * LP is processed at current node -> we can add violated bitvar constraints to the LP */

   *result = SCIP_FEASIBLE;

   /* step 1: check all useful bitvar constraints for feasibility */
   for( c = 0; c < nusefulconss; ++c )
   {
      debugMessage("LP enforcing bitvar constraint <%s>\n", SCIPconsGetName(conss[c]));
      CHECK_OKAY( separateCons(scip, conss[c], result) );
   }
   if( *result != SCIP_FEASIBLE )
      return SCIP_OKAY;

   /* step 2: check all obsolete bitvar constraints for feasibility */
   for( c = nusefulconss; c < nconss && *result == SCIP_FEASIBLE; ++c )
   {
      debugMessage("LP enforcing bitvar constraint <%s>\n", SCIPconsGetName(conss[c]));
      CHECK_OKAY( separateCons(scip, conss[c], result) );
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
DECL_CONSENFOPS(consEnfopsBitvar)
{
   Bool violated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* if the solution is infeasible anyway due to objective value, skip the enforcement */
   if( objinfeasible )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* check all bitvar constraints for feasibility */
   violated = FALSE;
   for( c = 0; c < nconss && !violated; ++c )
   {
      CHECK_OKAY( checkCons(scip, conss[c], NULL, TRUE, &violated) );
   }

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
DECL_CONSCHECK(consCheckBitvar)
{
   Bool violated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* check all bitvar constraints for feasibility */
   violated = FALSE;
   for( c = 0; c < nconss && !violated; ++c )
   {
      CHECK_OKAY( checkCons(scip, conss[c], sol, checklprows, &violated) );
   }

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
DECL_CONSPROP(consPropBitvar)
{
   Bool infeasible;
   int nfixedvars;
   int nchgbds;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* propagate all useful bitvar constraints */
   infeasible = FALSE;
   nfixedvars = 0;
   nchgbds = 0;
   for( c = 0; c < nusefulconss && !infeasible; ++c )
   {
      CHECK_OKAY( propagateCons(scip, conss[c], &nfixedvars, &nchgbds, NULL, &infeasible) );
   }

   /* adjust the result */
   if( infeasible )
      *result = SCIP_CUTOFF;
   else if( nfixedvars > 0 || nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
DECL_CONSPRESOL(consPresolBitvar)
{
   Bool infeasible;
   int nactfixedvars;
   int nactchgbds;
   int nactdelconss;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* propagate all bitvar constraints */
   nactfixedvars = 0;
   nactchgbds = 0;
   nactdelconss = 0;
   infeasible = FALSE;
   for( c = 0; c < nconss && !infeasible; ++c )
   {
      CHECK_OKAY( propagateCons(scip, conss[c], &nactfixedvars, &nactchgbds, &nactdelconss, &infeasible) );
   }

   /* adjust the result */
   if( infeasible )
      *result = SCIP_CUTOFF;
   else if( nactfixedvars > 0 || nactchgbds > 0 || nactdelconss > 0 )
   {
      *result = SCIP_SUCCESS;
      (*nfixedvars) += nactfixedvars;
      (*nchgbds) += nactchgbds;
      (*ndelconss) += nactdelconss;
   }

   return SCIP_OKAY;
}


/** conflict variable resolving method of constraint handler */
#define consRescvarBitvar NULL


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockBitvar)
{
   CONSDATA* consdata;
   int i;

   /* every change in variable values alters the feasibility of the bitvar equations:
    * we have to lock rounding in both directions for positive and negative constraint
    */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->bits != NULL);
   assert(consdata->words != NULL);

   /* lock bit variables */
   for( i = 0; i < consdata->nbits; ++i )
   {
      SCIPvarLock(consdata->bits[i], nlockspos + nlocksneg, nlockspos + nlocksneg);
   }

   /* lock word variables */
   for( i = 0; i < consdata->nwords; ++i )
   {
      SCIPvarLock(consdata->words[i], nlockspos + nlocksneg, nlockspos + nlocksneg);
   }

   return SCIP_OKAY;
}


/** variable rounding unlock method of constraint handler */
static
DECL_CONSUNLOCK(consUnlockBitvar)
{
   CONSDATA* consdata;
   int i;

   /* every change in variable values alters the feasibility of the bitvar equations:
    * we have to unlock rounding in both directions for positive and negative constraint
    */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->bits != NULL);
   assert(consdata->words != NULL);

   /* unlock bit variables */
   for( i = 0; i < consdata->nbits; ++i )
   {
      SCIPvarUnlock(consdata->bits[i], nunlockspos + nunlocksneg, nunlockspos + nunlocksneg);
   }

   /* unlock word variables */
   for( i = 0; i < consdata->nwords; ++i )
   {
      SCIPvarUnlock(consdata->words[i], nunlockspos + nunlocksneg, nunlockspos + nunlocksneg);
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#define consActiveBitvar NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveBitvar NULL


/** constraint enabling notification method of constraint handler */
#define consEnableBitvar NULL


/** constraint disabling notification method of constraint handler */
#define consDisableBitvar NULL




/*
 * bitvar event handler methods
 */

static
DECL_EVENTEXEC(eventExecBitvar)
{
   CONSDATA* consdata;

   consdata = (CONSDATA*)eventdata;
   assert(consdata != NULL);

   consdata->propagated = FALSE;

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for bitvar constraints and includes it in SCIP */
RETCODE SCIPincludeConsHdlrBitvar(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create event handler for bound change events */
   CHECK_OKAY( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
                  NULL, NULL, NULL,
                  NULL, eventExecBitvar,
                  NULL) );

   /* create bitvar constraint handler data */
   CHECK_OKAY( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConsHdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_NEEDSCONS,
                  consFreeBitvar, consInitBitvar, consExitBitvar, consSolstartBitvar,
                  consDeleteBitvar, consTransBitvar, consInitlpBitvar,
                  consSepaBitvar, consEnfolpBitvar, consEnfopsBitvar, consCheckBitvar, 
                  consPropBitvar, consPresolBitvar, consRescvarBitvar,
                  consLockBitvar, consUnlockBitvar,
                  consActiveBitvar, consDeactiveBitvar, 
                  consEnableBitvar, consDisableBitvar,
                  conshdlrdata) );

   return SCIP_OKAY;
}

/** creates and captures a bitvar constraint
 *  Warning! Either the bitvar should be short, or the objective value should be zero, because the objective
 *  value of the most significant bit in the variable would be 2^(nbits-1)*obj
 */
RETCODE SCIPcreateConsBitvar(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nbits,              /**< number of bits in the bitvar */
   Real             obj,                /**< objective value of bitvar variable */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             removeable          /**< should the constraint be removed from the LP due to aging or cleanup? */
   )
{
   CONSHDLR* conshdlr;
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;

   const Bool check = TRUE;       /* bitvar constraints must always be checked for feasibility */
   const Bool local = FALSE;      /* bitvar constraints are never local, because they represent problem variables */
   const Bool modifiable = FALSE; /* bitvar constraints are never modifiable, because they represent problem variables */

   /* find the bitvar constraint handler */
   conshdlr = SCIPfindConsHdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("bitvar constraint handler not found");
      return SCIP_PLUGINNOTFOUND;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create constraint data */
   CHECK_OKAY( consdataCreate(scip, &consdata, nbits) );
   CHECK_OKAY( consdataCreateVars(scip, consdata, conshdlrdata->eventhdlr, name, obj) );

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, 
                  initial, separate, enforce, check, propagate, local, modifiable, removeable) );

   return SCIP_OKAY;
}

/** creates and captures a constant bitvar constraint with constant given as bit vector
 *  Warning! Either the bitvar should be short, or the objective value should be zero, because the objective
 *  value of the most significant bit in the variable would be 2^(nbits-1)*obj
 */
RETCODE SCIPcreateConsBitconst(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nbits,              /**< number of bits in the bitvar */
   Real             obj,                /**< objective value of bitvar constant */
   Bool*            fixedbits           /**< array with the fixed bit values of constant (least significant first) */
   )
{
   CONSHDLR* conshdlr;
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;
   Bool infeasible;
   int b;

   const Bool initial = FALSE;     /* constant bitvar constraint doesn't have to be checked, propagated, ... */
   const Bool separate = FALSE;
   const Bool enforce = FALSE;
   const Bool propagate = FALSE;
   const Bool removeable = FALSE;
   const Bool check = FALSE;
   const Bool local = FALSE;
   const Bool modifiable = FALSE;

   assert(fixedbits != NULL);

   /* find the bitvar constraint handler */
   conshdlr = SCIPfindConsHdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("bitvar constraint handler not found");
      return SCIP_PLUGINNOTFOUND;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create constraint data */
   CHECK_OKAY( consdataCreate(scip, &consdata, nbits) );
   consdata->constant = TRUE;
   CHECK_OKAY( consdataCreateVars(scip, consdata, conshdlrdata->eventhdlr, name, obj) );

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, 
                  initial, separate, enforce, check, propagate, local, modifiable, removeable) );

   /* fix all variables according to the given bits */
   for( b = 0; b < consdata->nbits; ++b )
   {
      CHECK_OKAY( SCIPfixVar(scip, consdata->bits[b], (Real)(fixedbits[b]), &infeasible) );
      assert(!infeasible);
   }

   /* use propagation to fix the word to the corresponding value */
   CHECK_OKAY( propagateCons(scip, *cons, NULL, NULL, NULL, &infeasible) );
   assert(!infeasible);

   return SCIP_OKAY;
}

/** creates and captures a constant bitvar constraint with constant parsed from a string
 *  Warning! Either the bitvar should be short, or the objective value should be zero, because the objective
 *  value of the most significant bit in the variable would be 2^(nbits-1)*obj
 */
RETCODE SCIPcreateConsBitconstString(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nbits,              /**< number of bits in the bitvar */
   Real             obj,                /**< objective value of bitvar constant */
   const char*      cstring             /**< constant given as a string */
   )
{
   Bool* fixedbits;
   char msg[MAXSTRLEN];

   assert(cstring != NULL);

   /* get memory for constant as bit array */
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &fixedbits, nbits) );

   /* parse the string in a bit array */
   switch( *cstring )
   {
   case 'b':
   case 'B':
      CHECK_OKAY( parseBitString(cstring+1, nbits, fixedbits) );
      break;

   default:
      sprintf(msg, "invalid base character <%c> in bit constant string", *cstring);
      errorMessage(msg);
      return SCIP_PARSEERROR;
   }

   /* create the bitvar constant */
   CHECK_OKAY( SCIPcreateConsBitconst(scip, cons, name, nbits, obj, fixedbits) );

   /* release memory buffer */
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &fixedbits) );

   return SCIP_OKAY;
}

/** gets number of bits in bitvar */
int SCIPgetNBitsBitvar(
   CONS*            cons                /**< bitvar constraint */
   )
{
   CONSDATA* consdata;

   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nbits;
}

/** gets array with bits of bitvar, sorted least significant bit first */
VAR** SCIPgetBitsBitvar(
   CONS*            cons                /**< bitvar constraint */
   )
{
   CONSDATA* consdata;

   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->bits;
}

/** gets variable for single bit in bitvar */
VAR* SCIPgetBitBitvar(
   CONS*            cons,               /**< bitvar constraint */
   int              bit                 /**< bit number to get variable for */
   )
{
   CONSDATA* consdata;

   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= bit && bit < consdata->nbits);

   return consdata->bits[bit];
}

/** gets number of words in bitvar */
int SCIPgetNWordsBitvar(
   CONS*            cons                /**< bitvar constraint */
   )
{
   CONSDATA* consdata;

   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nwords;
}

/** gets array with words of bitvar, sorted least significant word first */
VAR** SCIPgetWordsBitvar(
   CONS*            cons                /**< bitvar constraint */
   )
{
   CONSDATA* consdata;

   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->words;
}

/** gets variable for single word in bitvar */
VAR* SCIPgetWordBitvar(
   CONS*            cons,               /**< bitvar constraint */
   int              word                /**< word number to get variable for */
   )
{
   CONSDATA* consdata;

   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= word && word < consdata->nwords);

   return consdata->words[word];
}

/** gets number of bits in a given word of a bitvar */
int SCIPgetNWordBitsBitvar(
   CONS*            cons,               /**< bitvar constraint */
   int              word                /**< word number */
   )
{
   CONSDATA* consdata;

   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(word >= 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( word >= consdata->nwords )
      return 0;
   else
      return wordSize(consdata, word);
}

/** returns the number of bits of the given word */
int SCIPgetWordSizeBitvar(
   CONS*            cons,               /**< bitvar constraint */
   int              word                /**< word number */
   )
{
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   return wordSize(SCIPconsGetData(cons), word);
}

/** returns the number of different values the given word can store (2^#bits) */
int SCIPgetWordPowerBitvar(
   CONS*            cons,               /**< bitvar constraint */
   int              word                /**< word number */
   )
{
   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   return wordPower(SCIPconsGetData(cons), word);
}

/** returns whether the bitvar was created as a constant */
Bool SCIPisConstBitvar(
   CONS*            cons                /**< bitvar constraint */
   )
{
   CONSDATA* consdata;

   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->constant;
}

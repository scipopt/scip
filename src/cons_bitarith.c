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

/**@file   cons_bitarith.c
 * @brief  constraint handler for bitarith constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_bitarith.h"
#include "cons_bitvar.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "bitarith"
#define CONSHDLR_DESC          "arithmetic constraints on fixed-length bit string variables"
#define CONSHDLR_SEPAPRIORITY  +1900000
#define CONSHDLR_ENFOPRIORITY  - 510000
#define CONSHDLR_CHECKPRIORITY - 510000
#define CONSHDLR_SEPAFREQ             1
#define CONSHDLR_PROPFREQ             1
#define CONSHDLR_NEEDSCONS         TRUE

#define EVENTHDLR_NAME         "bitarith"
#define EVENTHDLR_DESC         "bound change event handler for bitarith constraints"

#define BITVAR_NAME            "bitvar"




/*
 * Data structures
 */

/** arity of each arithmetic bit operation */
static const
int bitarithArity[SCIP_NBITARITHTYPES] =
{
   /*SCIP_BITARITHTYPE_ADD*/ 2,
   /*SCIP_BITARITHTYPE_SHL*/ 2,
   /*SCIP_BITARITHTYPE_EQ*/  2,
   /*SCIP_BITARITHTYPE_NOT*/ 1
};

/** constraint data for bitarith constraints */
struct ConsData
{
   CONS*            operand1;           /**< first (left) operand in operation (x) */
   CONS*            operand2;           /**< second (right) operand in operation (y) */
   CONS*            resultant;          /**< result of operation (z) */
   VAR**            vars;               /**< internal variables */
   ROW**            rows;               /**< LP rows representing linear relaxation of constraint */
   int              nvars;              /**< number of internal variables */
   int              nrows;              /**< number of LP rows used for representing linear relaxation of constraints */
   unsigned int     bitarithtype:2;     /**< type of arithmetic bit operation */
   unsigned int     propagated:1;       /**< is constraint already preprocessed/propagated? */
};

/** constraint handler data */
struct ConsHdlrData
{
   EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
};

/** state of a single bit variable: fixed to zero, unfixed, or fixed to one */
enum Bitstate
{
   FIXEDZERO = -1,                      /**< variable is fixed to zero; value has to be -1 to work with propagation! */
   UNFIXED   =  0,                      /**< variable is not fixed; value has to be 0 to work with propagation! */
   FIXEDONE  = +1                       /**< variable is fixed to one; value has to be +1 to work with propagation! */
};
typedef enum Bitstate BITSTATE;




/*
 * Local methods
 */

/** validates that the operands and resultant bitvars fit with the given arithmetic bit operation */
static
RETCODE validateArith(
   BITARITHTYPE     bitarithtype,       /**< type of arithmetic bit operation */
   CONS*            operand1,           /**< bitvar constraint: first (left) operand in operation (x) */
   CONS*            operand2,           /**< bitvar constraint: second (right) operand in operation (y) */
   CONS*            resultant           /**< bitvar constraint: result of operation (z) */
   )
{
   /* check the existence of the arithmetic operation */
   if( bitarithtype >= SCIP_NBITARITHTYPES )
   {
      errorMessage("unknown arithmetic bit operation");
      return SCIP_INVALIDDATA;
   }

   /* check the arity of the arithmetic operation */
   if( (bitarithArity[bitarithtype] < 1) != (operand1 == NULL)
      || (bitarithArity[bitarithtype] < 2) != (operand2 == NULL)
      || resultant == NULL )
   {
      errorMessage("invalid number of operands given for bitarith constraint");
      return SCIP_INVALIDDATA;
   }

   /* check if the operand and resultant constraints are already added to the problem */
   if( !SCIPconsIsInProb(operand1) && !SCIPisConstBitvar(operand1) )
   {
      errorMessage("first operand of bitarith constraint must either be constant or been added to the problem");
      return SCIP_INVALIDDATA;
   }
   if( !SCIPconsIsInProb(operand2) && !SCIPisConstBitvar(operand2) )
   {
      errorMessage("second operand of bitarith constraint must either be constant or been added to the problem");
      return SCIP_INVALIDDATA;
   }
   if( !SCIPconsIsInProb(resultant) && !SCIPisConstBitvar(resultant) )
   {
      errorMessage("resultant of bitarith constraint must either be constant or been added to the problem");
      return SCIP_INVALIDDATA;
   }

   /* check the lengths of the operands and resultant */
   switch( bitarithtype )
   {
   case SCIP_BITARITHTYPE_ADD:
      if( SCIPgetNBitsBitvar(operand1) > SCIPgetNBitsBitvar(resultant)
         || SCIPgetNBitsBitvar(operand2) > SCIPgetNBitsBitvar(resultant) )
      {
         errorMessage("operand size must not be larger than resultant size in bit-add constraint");
         return SCIP_INVALIDDATA;
      }
      break;

   case SCIP_BITARITHTYPE_SHL:
   case SCIP_BITARITHTYPE_EQ:
   case SCIP_BITARITHTYPE_NOT:
      errorMessage("not implemented yet");
      abort();
      break;

   default:
      errorMessage("invalid bit arithmetic type");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** creates constaint handler data for bitarith constraint handler */
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
      errorMessage("event handler for bitarith constraints not found");
      return SCIP_PLUGINNOTFOUND;
   }

   return SCIP_OKAY;
}

/** frees constaint handler data for bitarith constraint handler */
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

/** catches all events associated with the given bitarith constraint data */
static
RETCODE consdataCatchEvents(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< bitarith constraint data */
   EVENTHDLR*       eventhdlr           /**< event handler for bound change events */
   )
{
   VAR** bits;
   int nbits;
   int b;
   int v;

   /* catch bound tighten events on all bit variables in operands and resultant */
   if( consdata->operand1 != NULL )
   {
      bits = SCIPgetBitsBitvar(consdata->operand1);
      nbits = SCIPgetNBitsBitvar(consdata->operand1);
      for( b = 0; b < nbits; ++b )
      {
         CHECK_OKAY( SCIPcatchVarEvent(scip, bits[b], SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr, (EVENTDATA*)consdata) );
      }
   }
   if( consdata->operand2 != NULL )
   {
      bits = SCIPgetBitsBitvar(consdata->operand2);
      nbits = SCIPgetNBitsBitvar(consdata->operand2);
      for( b = 0; b < nbits; ++b )
      {
         CHECK_OKAY( SCIPcatchVarEvent(scip, bits[b], SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr, (EVENTDATA*)consdata) );
      }
   }
   if( consdata->resultant != NULL )
   {
      bits = SCIPgetBitsBitvar(consdata->resultant);
      nbits = SCIPgetNBitsBitvar(consdata->resultant);
      for( b = 0; b < nbits; ++b )
      {
         CHECK_OKAY( SCIPcatchVarEvent(scip, bits[b], SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr, (EVENTDATA*)consdata) );
      }
   }

   /* catch bound tighten events on internal variables */
   for( v = 0; v < consdata->nvars; ++v )
   {
      CHECK_OKAY( SCIPcatchVarEvent(scip, consdata->vars[v], SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr,
                     (EVENTDATA*)consdata) );
   }

   return SCIP_OKAY;
}

/** drops all events associated with the given bitarith constraint data */
static
RETCODE consdataDropEvents(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< bitarith constraint data */
   EVENTHDLR*       eventhdlr           /**< event handler for bound change events */
   )
{
   VAR** bits;
   int nbits;
   int b;
   int v;

   /* drop bound tighten events on all bit variables in operands and resultant */
   if( consdata->operand1 != NULL )
   {
      bits = SCIPgetBitsBitvar(consdata->operand1);
      nbits = SCIPgetNBitsBitvar(consdata->operand1);
      for( b = 0; b < nbits; ++b )
      {
         CHECK_OKAY( SCIPdropVarEvent(scip, bits[b], eventhdlr, (EVENTDATA*)consdata) );
      }
   }
   if( consdata->operand2 != NULL )
   {
      bits = SCIPgetBitsBitvar(consdata->operand2);
      nbits = SCIPgetNBitsBitvar(consdata->operand2);
      for( b = 0; b < nbits; ++b )
      {
         CHECK_OKAY( SCIPdropVarEvent(scip, bits[b], eventhdlr, (EVENTDATA*)consdata) );
      }
   }
   if( consdata->resultant != NULL )
   {
      bits = SCIPgetBitsBitvar(consdata->resultant);
      nbits = SCIPgetNBitsBitvar(consdata->resultant);
      for( b = 0; b < nbits; ++b )
      {
         CHECK_OKAY( SCIPdropVarEvent(scip, bits[b], eventhdlr, (EVENTDATA*)consdata) );
      }
   }

   /* drop bound tighten events on internal variables */
   for( v = 0; v < consdata->nvars; ++v )
   {
      CHECK_OKAY( SCIPdropVarEvent(scip, consdata->vars[v], eventhdlr, (EVENTDATA*)consdata) );
   }

   return SCIP_OKAY;
}

/** creates a bitarith constraint data object and captures the given operand and resultant bitvar constraints */
static
RETCODE consdataCreate(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to store the bitarith constraint data */
   BITARITHTYPE     bitarithtype,       /**< type of arithmetic bit operation */
   CONS*            operand1,           /**< bitvar constraint: first (left) operand in operation (x) */
   CONS*            operand2,           /**< bitvar constraint: second (right) operand in operation (y) */
   CONS*            resultant           /**< bitvar constraint: result of operation (z) */
   )
{
   assert(consdata != NULL);
   assert(bitarithtype < SCIP_NBITARITHTYPES);
   assert((bitarithArity[bitarithtype] < 1) == (operand1 == NULL));
   assert((bitarithArity[bitarithtype] < 2) == (operand2 == NULL));
   assert(resultant != NULL);

   /* create constraint data */
   CHECK_OKAY( SCIPallocBlockMemory(scip, consdata) );
   (*consdata)->operand1 = operand1;
   (*consdata)->operand2 = operand2;
   (*consdata)->resultant = resultant;
   (*consdata)->vars = NULL;
   (*consdata)->rows = NULL;
   (*consdata)->nvars = 0;
   (*consdata)->nrows = 0;
   (*consdata)->bitarithtype = bitarithtype;
   (*consdata)->propagated = FALSE;

   /* capture operand and resultant constraints */
   if( (*consdata)->operand1 != NULL )
   {
      assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr((*consdata)->operand1)), BITVAR_NAME) == 0);
      CHECK_OKAY( SCIPcaptureCons(scip, (*consdata)->operand1) );
   }
   if( (*consdata)->operand2 != NULL )
   {
      assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr((*consdata)->operand2)), BITVAR_NAME) == 0);
      CHECK_OKAY( SCIPcaptureCons(scip, (*consdata)->operand2) );
   }
   if( (*consdata)->resultant != NULL )
   {
      assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr((*consdata)->resultant)), BITVAR_NAME) == 0);
      CHECK_OKAY( SCIPcaptureCons(scip, (*consdata)->resultant) );
   }

   return SCIP_OKAY;
}

/** frees a bitarith constraint data object and releases corresponding variables */
static
RETCODE consdataFree(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA**       consdata,           /**< pointer to the bitarith constraint data */
   EVENTHDLR*       eventhdlr           /**< event handler for bound change events */
   )
{
   int i;

   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* drop events on transformed constraints */
   if( SCIPisTransformed(scip) )
   {
      CHECK_OKAY( consdataDropEvents(scip, *consdata, eventhdlr) );
   }

   /* release operand and resultant constraints */
   if( (*consdata)->operand1 != NULL )
   {
      CHECK_OKAY( SCIPreleaseCons(scip, &(*consdata)->operand1) );
   }
   if( (*consdata)->operand2 != NULL )
   {
      CHECK_OKAY( SCIPreleaseCons(scip, &(*consdata)->operand2) );
   }
   if( (*consdata)->resultant != NULL )
   {
      CHECK_OKAY( SCIPreleaseCons(scip, &(*consdata)->resultant) );
   }

   /* release internal variables */
   for( i = 0; i < (*consdata)->nvars; ++i )
   {
      CHECK_OKAY( SCIPreleaseVar(scip, &(*consdata)->vars[i]) );
   }

   /* release rows */
   for( i = 0; i < (*consdata)->nrows; ++i )
   {
      CHECK_OKAY( SCIPreleaseRow(scip, &(*consdata)->rows[i]) );
   }

   /* free memory */
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vars, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->rows, (*consdata)->nrows);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** replaces the operands and resultant constraints with their transformed counterparts, if existing */
static
RETCODE consdataTransformOperands(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata            /**< bitarith constraint data */
   )
{
   CONS* cons;

   assert(SCIPstage(scip) > SCIP_STAGE_PROBLEM);

   /* transform bitvar constraints, if they are not yet transformed */
   if( consdata->operand1 != NULL && !SCIPconsIsTransformed(consdata->operand1) )
   {
      CHECK_OKAY( SCIPtransformCons(scip, consdata->operand1, &cons) );
      CHECK_OKAY( SCIPreleaseCons(scip, &consdata->operand1) );
      consdata->operand1 = cons;
   }
   if( consdata->operand2 != NULL && !SCIPconsIsTransformed(consdata->operand2) )
   {
      CHECK_OKAY( SCIPtransformCons(scip, consdata->operand2, &cons) );
      CHECK_OKAY( SCIPreleaseCons(scip, &consdata->operand2) );
      consdata->operand2 = cons;
   }
   if( consdata->resultant != NULL && !SCIPconsIsTransformed(consdata->resultant) )
   {
      CHECK_OKAY( SCIPtransformCons(scip, consdata->resultant, &cons) );
      CHECK_OKAY( SCIPreleaseCons(scip, &consdata->resultant) );
      consdata->resultant = cons;
   }

   return SCIP_OKAY;
}

/** creates and adds the given number of internal binary variables for the given bitarith constraint data */
static
RETCODE consdataCreateVars(
   SCIP*            scip,               /**< SCIP data structure */
   CONSDATA*        consdata,           /**< bitarith constraint data */
   int              nvars,              /**< number of variables to create */
   const char*      name                /**< name prefix for variables */
   )
{
   char varname[MAXSTRLEN];
   int v;

   assert(consdata != NULL);

   debugMessage("creating %d variables with name prefix <%s> for bitarith constraint\n", nvars, name);

   CHECK_OKAY( SCIPallocBlockMemoryArray(scip, &consdata->vars, nvars) );
   consdata->nvars = nvars;

   for( v = 0; v < nvars; ++v )
   {
      sprintf(varname, "%s_%d", name, v);
      CHECK_OKAY( SCIPcreateVar(scip, &consdata->vars[v], varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, TRUE) );
      CHECK_OKAY( SCIPaddVar(scip, consdata->vars[v]) );
   }

   return SCIP_OKAY;
}

/** creates and adds the necessary internal variables for the given bitarith constraint data */
static
RETCODE createVars(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< bitarith constraint */
   )
{
   CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   switch( consdata->bitarithtype )
   {
   case SCIP_BITARITHTYPE_ADD:
      /* for each word in the resultant, one overflow binary variable is needed */
      CHECK_OKAY( consdataCreateVars(scip, consdata, SCIPgetNWordsBitvar(consdata->resultant), SCIPconsGetName(cons)) );
      break;

   case SCIP_BITARITHTYPE_SHL:
   case SCIP_BITARITHTYPE_EQ:
   case SCIP_BITARITHTYPE_NOT:
      errorMessage("not implemented yet");
      abort();
      break;

   default:
      errorMessage("invalid bit arithmetic type");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** creates LP relaxation for add constraint */
static
RETCODE createRowsAdd(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< bitarith constraint */
   )
{
   CONSDATA* consdata;
   VAR** op1words;
   VAR** op2words;
   VAR** reswords;
   char rowname[MAXSTRLEN];
   int nop1words;
   int nop2words;
   int nreswords;
   int w;

   /* each word has one associated row: 
    *   resultant[w] = operand1[w] + operand2[w] + overflow[w-1] - wordpower[w] * overflow[w]
    */

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->bitarithtype == SCIP_BITARITHTYPE_ADD);
   assert(consdata->rows == NULL);
   assert(consdata->nrows == 0);

   /* get the word variables of the bitvars */
   op1words = SCIPgetWordsBitvar(consdata->operand1);
   nop1words = SCIPgetNWordsBitvar(consdata->operand1);
   assert(nop1words == 0 || op1words != NULL);
   op2words = SCIPgetWordsBitvar(consdata->operand2);
   nop2words = SCIPgetNWordsBitvar(consdata->operand2);
   assert(nop2words == 0 || op2words != NULL);
   reswords = SCIPgetWordsBitvar(consdata->resultant);
   nreswords = SCIPgetNWordsBitvar(consdata->resultant);
   assert(nreswords == 0 || reswords != NULL);
   assert(nreswords == consdata->nvars);

   /* get memory for rows */
   consdata->nrows = nreswords;
   CHECK_OKAY( SCIPallocBlockMemoryArray(scip, &consdata->rows, consdata->nrows) );
   
   /* create the rows:
    *   - resultant[w] + operand1[w] + operand2[w] + overflow[w-1] - wordpower[w] * overflow[w] == 0.0,
    *   with operand1[w] == 0.0 and operand2[w] == 0.0, if the word doesn't exist (because operator is too small),
    *   and overflow[-1] == 0.0
    */
   for( w = 0; w < nreswords; ++w )
   {
      sprintf(rowname, "%s_r%d", SCIPconsGetName(cons), w);
      CHECK_OKAY( SCIPcreateEmptyRow(scip, &consdata->rows[w], rowname, 0.0, 0.0,
                     SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );

      /* - resultant[w] */
      CHECK_OKAY( SCIPaddVarToRow(scip, consdata->rows[w], reswords[w], -1.0) );
      /* + operand1[w] */
      if( w < nop1words )
      {
         CHECK_OKAY( SCIPaddVarToRow(scip, consdata->rows[w], op1words[w], +1.0) );
      }
      /* + operand2[w] */
      if( w < nop2words )
      {
         CHECK_OKAY( SCIPaddVarToRow(scip, consdata->rows[w], op2words[w], +1.0) );
      }
      /* + overflow[w-1] */
      if( w != 0 )
      {
         CHECK_OKAY( SCIPaddVarToRow(scip, consdata->rows[w], consdata->vars[w-1], +1.0) );
      }
      /* - wordpower[w] * overflow[w] */
      CHECK_OKAY( SCIPaddVarToRow(scip, consdata->rows[w], consdata->vars[w], 
                     -SCIPgetWordPowerBitvar(consdata->resultant, w)) );

      debugMessage("created row <%s> for bitarith add constraint <%s> word %d\n", 
         SCIProwGetName(consdata->rows[w]), SCIPconsGetName(cons), w);
      debug(SCIPprintRow(scip, consdata->rows[w], NULL));
   }

   return SCIP_OKAY;
}

/** creates LP relaxation for bitarith constraint */
static
RETCODE createRows(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< bitarith constraint */
   )
{
   CONSDATA* consdata;
   int r;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->rows == NULL);
   assert(consdata->nrows == 0);

   /* create the LP relaxation, depending on the arithmetic operation */
   switch( consdata->bitarithtype )
   {
   case SCIP_BITARITHTYPE_ADD:
      CHECK_OKAY( createRowsAdd(scip, cons) );
      break;

   case SCIP_BITARITHTYPE_SHL:
   case SCIP_BITARITHTYPE_EQ:
   case SCIP_BITARITHTYPE_NOT:
      errorMessage("not implemented yet");
      abort();
      break;

   default:
      errorMessage("invalid bit arithmetic type");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** checks a single word of an add constraint for feasibility of given solution */
static
RETCODE checkAddWord(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitarith constraint */
   int              word,               /**< word number to check */
   SOL*             sol,                /**< solution to be checked, or NULL for actual solution */
   Bool             checklprows,        /**< has bitvar constraint to be checked, if it is already in current LP? */
   Real*            violation,          /**< pointer to store the constraint's violation, or NULL */
   Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   CONSDATA* consdata;
   Real sum;

   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == SCIPgetNWordsBitvar(consdata->resultant));
   assert(0 <= word && word < consdata->nvars);

   debugMessage("checking bitarith add constraint <%s> at word %d\n", SCIPconsGetName(cons), word);

   *violated = FALSE;

   if( !checklprows && consdata->rows != NULL && consdata->rows[word] != NULL && SCIProwIsInLP(consdata->rows[word]) )
      return SCIP_OKAY;

   /* we have to check the constraint
    *   resultant[w] == operand1[w] + operand2[w] + overflow[w-1] - wordpower[w] * overflow[w],
    *   with operand1[w] == 0.0 and operand2[w] == 0.0, if the word doesn't exist (because operator is too small),
    *   and overflow[-1] == 0.0
    */
   sum = -SCIPgetSolVal(scip, sol, SCIPgetWordBitvar(consdata->resultant, word));
   sum -= SCIPgetWordPowerBitvar(consdata->resultant, word) * SCIPgetSolVal(scip, sol, consdata->vars[word]);
   if( word < SCIPgetNWordsBitvar(consdata->operand1) )
      sum += SCIPgetSolVal(scip, sol, SCIPgetWordBitvar(consdata->operand1, word));
   if( word < SCIPgetNWordsBitvar(consdata->operand2) )
      sum += SCIPgetSolVal(scip, sol, SCIPgetWordBitvar(consdata->operand2, word));
   if( word > 0 )
      sum += SCIPgetSolVal(scip, sol, consdata->vars[word-1]);
   
   /* check violation: sum should be zero -> violation = ABS(sum) */
   if( violation != NULL )
      *violation = ABS(sum);
   *violated = !SCIPisFeasZero(scip, sum);

   return SCIP_OKAY;
}

/** checks an add constraint for feasibility of given solution */
static
RETCODE checkAdd(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitarith constraint */
   SOL*             sol,                /**< solution to be checked, or NULL for actual solution */
   Bool             checklprows,        /**< has bitvar constraint to be checked, if it is already in current LP? */
   Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   CONSDATA* consdata;
   int w;

   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == SCIPgetNWordsBitvar(consdata->resultant));

   /* check each word in the resultant bitvar */
   *violated = FALSE;
   for( w = 0; w < consdata->nvars && !(*violated); ++w )
   {
      CHECK_OKAY( checkAddWord(scip, cons, w, sol, checklprows, NULL, violated) );
   }

   return SCIP_OKAY;
}

/** checks a bitarith constraint for feasibility of given solution */
static
RETCODE checkCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitarith constraint */
   SOL*             sol,                /**< solution to be checked, or NULL for actual solution */
   Bool             checklprows,        /**< has bitvar constraint to be checked, if it is already in current LP? */
   Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   switch( consdata->bitarithtype )
   {
   case SCIP_BITARITHTYPE_ADD:
      CHECK_OKAY( checkAdd(scip, cons, sol, checklprows, violated) );
      break;

   case SCIP_BITARITHTYPE_SHL:
   case SCIP_BITARITHTYPE_EQ:
   case SCIP_BITARITHTYPE_NOT:
      errorMessage("not implemented yet");
      abort();
      break;

   default:
      errorMessage("invalid bit arithmetic type");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** adds single row of bitarith constraint as cut to the LP */
static
RETCODE addCut(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitarith constraint */
   int              rownum,             /**< row number to add as cut */
   Real             cutscore            /**< basic score value for the added cuts */
   )
{
   CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert((consdata->nrows == 0) == (consdata->rows == NULL));

   /* create the corresponding rows, if not yet existing */
   if( consdata->rows == NULL )
   {
      /* convert consdata object into LP row */
      CHECK_OKAY( createRows(scip, cons) );
   }
   assert(consdata->rows != NULL);
   assert(consdata->nrows != 0);
   assert(0 <= rownum && rownum < consdata->nrows);

   debugMessage("adding cut <%s> from bitarith constraint <%s> row %d\n",
      SCIProwGetName(consdata->rows[rownum]), SCIPconsGetName(cons), rownum);
   debug(SCIPprintRow(scip, consdata->rows[rownum], NULL));

   /* insert LP row as cuts */
   assert(!SCIProwIsInLP(consdata->rows[rownum]));
   CHECK_OKAY( SCIPaddCut(scip, consdata->rows[rownum], 
                  cutscore/SCIProwGetNorm(consdata->rows[rownum])/(SCIProwGetNNonz(consdata->rows[rownum])+1)) );
   
   return SCIP_OKAY;
}

/** adds all non-LP rows of bitarith constraint as cuts to the LP */
static
RETCODE addAllCuts(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitarith constraint */
   Real             cutscore            /**< basic score value for the added cuts */
   )
{
   CONSDATA* consdata;
   int r;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* create the corresponding rows, if not yet existing */
   if( consdata->rows == NULL )
   {
      /* convert consdata object into LP row */
      CHECK_OKAY( createRows(scip, cons) );
   }
   assert(consdata->rows != NULL);
   assert(consdata->nrows != 0);

   /** add all rows to the LP */
   for( r = 0; r < consdata->nrows; ++r )
   {
      if( !SCIProwIsInLP(consdata->rows[r]) )
      {
         CHECK_OKAY( addCut(scip, cons, r, cutscore) );
      }
   }

   return SCIP_OKAY;
}

/** separates add constraint */
static
RETCODE separateAdd(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitarith constraint */
   RESULT*          result              /**< pointer to store result of separation */
   )
{
   CONSDATA* consdata;
   Real violation;
   Bool violated;
   int w;

   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == SCIPgetNWordsBitvar(consdata->resultant));

   for( w = 0; w < consdata->nvars; ++w )
   {
      CHECK_OKAY( checkAddWord(scip, cons, w, NULL, FALSE, &violation, &violated) );
      
      if( violated )
      {
         /* insert LP row as cut */
         CHECK_OKAY( addCut(scip, cons, w, violation) );
         *result = SCIP_SEPARATED;
      }
   }

   return SCIP_OKAY;
}

/** separates bitarith constraint */
static
RETCODE separateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitarith constraint */
   RESULT*          result              /**< pointer to store result of separation */
   )
{
   CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   switch( consdata->bitarithtype )
   {
   case SCIP_BITARITHTYPE_ADD:
      CHECK_OKAY( separateAdd(scip, cons, result) );
      break;

   case SCIP_BITARITHTYPE_SHL:
   case SCIP_BITARITHTYPE_EQ:
   case SCIP_BITARITHTYPE_NOT:
      errorMessage("not implemented yet");
      abort();
      break;

   default:
      errorMessage("invalid bit arithmetic type");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** sets value in bitstate to given state, and sets *infeasible = TRUE, if bitstate is already fixed to opposite value;
 *  returns TRUE iff bitstate was unfixed before
 */
static
Bool fixState(
   BITSTATE*        bitstate,           /**< bitstate array */
   int              bit,                /**< bit number to set the new state for */
   BITSTATE         newstate,           /**< new state for bit */
   Bool*            infeasible          /**< pointer to store TRUE, if new state is infeasible */
   )
{
   assert(bitstate != NULL);
   assert(newstate == FIXEDZERO || newstate == FIXEDONE);
   assert(infeasible != NULL);

   if( bitstate[bit] == UNFIXED )
   {
      bitstate[bit] = newstate;
      return TRUE;
   }
   else
   {
      *infeasible = (bitstate[bit] != newstate);
      return FALSE;
   }
}

/** applies (and counts) a given binary variable fixing */
static
RETCODE deduceBinVar(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitarith constraint */
   VAR*             var,                /**< variable to deduce value for */
   Real             val,                /**< deduced value */
   int*             nfixedvars,         /**< pointer to add up the number of fixed variables, or NULL */
   Bool*            infeasible          /**< pointer to store TRUE, if fixed value is infeasible */
   )
{
   assert(infeasible != NULL);

   if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
   {
      /* variable already fixed */
      if( !SCIPisEQ(scip, SCIPvarGetLbLocal(var), val) )
         *infeasible = TRUE;
   }
   else
   {
      debugMessage("bitarith constraint <%s>: fixing variable <%s> to %g\n", 
         SCIPconsGetName(cons), SCIPvarGetName(var), val);
      
      CHECK_OKAY( SCIPinferBinVar(scip, var, val, cons) );
      if( nfixedvars != NULL )
         (*nfixedvars)++;
   }

   return SCIP_OKAY;
}

/** aggregates a binary variable $x$ to be equal to a binary variable $y$ */
static
RETCODE aggregateBinVarsPos(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitarith constraint */
   VAR*             var,                /**< variable x in aggregation $x - y = 0$ */
   VAR*             aggvar,             /**< variable y in aggregation $x - y = 0$ */
   int*             naggrvars,          /**< pointer to add up the number of aggregated variables */
   Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   Bool*            redundant           /**< pointer to store whether the equality is (now) redundant */
   )
{
   Bool aggregated;

   assert(naggrvars != NULL);

   debugMessage("bitarith constraint <%s>: aggregating variables <%s> - <%s> = 0\n", 
      SCIPconsGetName(cons), SCIPvarGetName(var), SCIPvarGetName(aggvar));

   CHECK_OKAY( SCIPaggregateVars(scip, var, aggvar, 1.0, -1.0, 0.0, infeasible, redundant, &aggregated) );
   if( aggregated )
      (*naggrvars)++;

   return SCIP_OKAY;
}

/** aggregates a binary variable $x$ to be equal to the negation of a binary variable $y$ */
static
RETCODE aggregateBinVarsNeg(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitarith constraint */
   VAR*             var,                /**< variable x in aggregation $x + y = 1$ */
   VAR*             aggvar,             /**< variable y in aggregation $x + y = 1$ */
   int*             naggrvars,          /**< pointer to add up the number of aggregated variables */
   Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   Bool*            redundant           /**< pointer to store whether the equality is (now) redundant */
   )
{
   Bool aggregated;

   assert(naggrvars != NULL);

   debugMessage("bitarith constraint <%s>: aggregating variables <%s> + <%s> = 1\n", 
      SCIPconsGetName(cons), SCIPvarGetName(var), SCIPvarGetName(aggvar));

   CHECK_OKAY( SCIPaggregateVars(scip, var, aggvar, 1.0, 1.0, 1.0, infeasible, redundant, &aggregated) );
   if( aggregated )
      (*naggrvars)++;

   return SCIP_OKAY;
}

/** aggregates variables of an add constraint with the rule:
 *    (1) if 2 of x_b, y_b, o_b, and z_b are fixed, one of the other two variables and o_(b+1) can be aggregated
 */
static
RETCODE aggregateAdd(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitarith constraint */
   BITSTATE*        op1state,           /**< current states of operand1 bits */
   BITSTATE*        op2state,           /**< current states of operand2 bits */
   BITSTATE*        ovrstate,           /**< current states of overflow bits */
   BITSTATE*        resstate,           /**< current states of resultant bits */
   VAR**            op1vars,            /**< binary variables representing operand1 */
   VAR**            op2vars,            /**< binary variables representing operand2 */
   VAR**            resvars,            /**< binary variables representing resultant */
   int*             naggrvars,          /**< pointer to add up the number of aggregated variables */
   Bool*            infeasible          /**< pointer to store whether the aggregation is infeasible */
   )
{
   CONSDATA* consdata;
   Bool op1fixed;
   Bool op2fixed;
   Bool ovrfixed;
   Bool resfixed;
   Bool redundant;
   int nresbits;
   int word;
   int wordbitstart;
   int nextwordbitstart;
   int b;

   assert(SCIPstage(scip) == SCIP_STAGE_PRESOLVING);
   assert(op1state != NULL);
   assert(op2state != NULL);
   assert(ovrstate != NULL);
   assert(resstate != NULL);
   assert(op1vars != NULL);
   assert(op2vars != NULL);
   assert(resvars != NULL);
   assert(naggrvars != NULL);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* get number of bits */
   nresbits = SCIPgetNBitsBitvar(consdata->resultant);

   /* if in the equation  x_b + y_b - o_b == z_b + 2*o_(b+1)  two of the variables x_b, y_b, o_b, and z_b are fixed,
    * the other ones and o_(b+1) can be aggregated:
    *  (1) x_b, y_b fixed:
    *    (a) x_b == y_b   =>    o_b == z_b        o_(b+1) fixed (done in propagation)
    *    (b) x_b != y_b   =>    o_b == 1 - z_b    o_(b+1) == 1 - z_b
    *  (2) x_b, o_b fixed:
    *    (a) x_b == o_b   =>    y_b == z_b        o_(b+1) fixed (done in propagation)
    *    (b) x_b != o_b   =>    y_b == 1 - z_b    o_(b+1) == 1 - z_b
    *  (3) x_b, z_b fixed:
    *    (a) x_b == z_b   =>    o_b == y_b        o_(b+1) == y_b
    *    (b) x_b != z_b   =>    o_b == 1 - y_b    o_(b+1) fixed (done in propagation)
    *  (4) y_b, o_b fixed:
    *    (a) y_b == o_b   =>    x_b == z_b        o_(b+1) fixed (done in propagation)
    *    (b) y_b != o_b   =>    x_b == 1 - z_b    o_(b+1) == 1 - z_b
    *  (5) y_b, z_b fixed:
    *    (a) y_b == z_b   =>    o_b == x_b        o_(b+1) == x_b
    *    (b) y_b != z_b   =>    o_b == 1 - x_b    o_(b+1) fixed (done in propagation)
    *  (6) o_b, z_b fixed:
    *    (a) o_b == z_b   =>    x_b == y_b        o_(b+1) == y_b
    *    (b) o_b != z_b   =>    x_b == 1 - y_b    o_(b+1) fixed (done in propagation)
    */
   word = 0;
   wordbitstart = 0;
   nextwordbitstart = SCIPgetWordSizeBitvar(consdata->resultant, 0);
   for( b = 0; b < nresbits && !(*infeasible); ++b )
   {
      /* check, if we are in the next word */
      if( b == nextwordbitstart )
      {
         word++;
         wordbitstart = nextwordbitstart;
         nextwordbitstart += SCIPgetWordSizeBitvar(consdata->resultant, word);
      }
      assert(wordbitstart < nextwordbitstart);
      assert(nextwordbitstart <= nresbits);
      assert(wordbitstart <= b && b < nextwordbitstart);

      /* check for possible aggregations */
      op1fixed = (op1state[b] != UNFIXED);
      op2fixed = (op2state[b] != UNFIXED);
      ovrfixed = (ovrstate[b] != UNFIXED);
      resfixed = (resstate[b] != UNFIXED);
      if( op1fixed + op2fixed + ovrfixed + resfixed == 2 )
      {
         if( op1fixed && op2fixed )
         {
            /* (1) x_b, y_b fixed */
            if( op1state[b] == op2state[b] )
            {
               /* (a) x_b == y_b   =>    o_b == z_b        o_(b+1) fixed (done in propagation) */
               assert(ovrstate[b+1] == op1state[b]);
               if( b == wordbitstart )
               {
                  CHECK_OKAY( aggregateBinVarsPos(scip, cons, consdata->vars[word], resvars[b], 
                                 naggrvars, infeasible, &redundant) );
               }
            }
            else
            {
               /* (b) x_b != y_b   =>    o_b == 1 - z_b    o_(b+1) == 1 - z_b */
               if( b == wordbitstart )
               {
                  CHECK_OKAY( aggregateBinVarsNeg(scip, cons, consdata->vars[word], resvars[b],
                                 naggrvars, infeasible, &redundant) );
               }
               if( b+1 == nextwordbitstart && !(*infeasible) )
               {
                  CHECK_OKAY( aggregateBinVarsNeg(scip, cons, consdata->vars[word+1], resvars[b],
                                 naggrvars, infeasible, &redundant) );
               }
            }
         }
         else if( op1fixed && ovrfixed )
         {
            /* (2) x_b, o_b fixed */
            if( op1state[b] == ovrstate[b] )
            {
               /* (a) x_b == o_b   =>    y_b == z_b        o_(b+1) fixed (done in propagation) */
               assert(ovrstate[b+1] == op1state[b]);
               CHECK_OKAY( aggregateBinVarsPos(scip, cons, op2vars[b], resvars[b], naggrvars, infeasible, &redundant) );
            }
            else
            {
               /* (b) x_b != o_b   =>    y_b == 1 - z_b    o_(b+1) == 1 - z_b */
               CHECK_OKAY( aggregateBinVarsNeg(scip, cons, op2vars[b], resvars[b], naggrvars, infeasible, &redundant) );
               if( b+1 == nextwordbitstart && !(*infeasible) )
               {
                  CHECK_OKAY( aggregateBinVarsNeg(scip, cons, consdata->vars[word+1], resvars[b],
                                 naggrvars, infeasible, &redundant) );
               }
            }
         }
         else if( op1fixed && resfixed )
         {
            /* (3) x_b, z_b fixed */
            if( op1state[b] == resstate[b] )
            {
               /* (a) x_b == z_b   =>    o_b == y_b        o_(b+1) == y_b */
               if( b == wordbitstart )
               {
                  CHECK_OKAY( aggregateBinVarsPos(scip, cons, consdata->vars[word], op2vars[b],
                                 naggrvars, infeasible, &redundant) );
               }
               if( b+1 == nextwordbitstart && !(*infeasible) )
               {
                  CHECK_OKAY( aggregateBinVarsPos(scip, cons, consdata->vars[word+1], op2vars[b],
                                 naggrvars, infeasible, &redundant) );
               }
            }
            else
            {
               /* (b) x_b != z_b   =>    o_b == 1 - y_b    o_(b+1) fixed (done in propagation) */
               assert(ovrstate[b+1] == op1state[b]);
               if( b == wordbitstart )
               {
                  CHECK_OKAY( aggregateBinVarsNeg(scip, cons, consdata->vars[word], op2vars[b],
                                 naggrvars, infeasible, &redundant) );
               }
            }
         }
         else if( op2fixed && ovrfixed )
         {
            /* (4) y_b, o_b fixed */
            if( op2state[b] == ovrstate[b] )
            {
               /* (a) y_b == o_b   =>    x_b == z_b        o_(b+1) fixed (done in propagation) */
               assert(ovrstate[b+1] == op2state[b]);
               CHECK_OKAY( aggregateBinVarsPos(scip, cons, op1vars[b], resvars[b], naggrvars, infeasible, &redundant) );
            }
            else
            {
               /* (b) y_b != o_b   =>    x_b == 1 - z_b    o_(b+1) == 1 - z_b */
               CHECK_OKAY( aggregateBinVarsNeg(scip, cons, op1vars[b], resvars[b], naggrvars, infeasible, &redundant) );
               if( b+1 == nextwordbitstart && !(*infeasible) )
               {
                  CHECK_OKAY( aggregateBinVarsNeg(scip, cons, consdata->vars[word+1], resvars[b],
                                 naggrvars, infeasible, &redundant) );
               }
            }
         }
         else if( op2fixed && resfixed )
         {
            /* (5) y_b, z_b fixed */
            if( op2state[b] == resstate[b] )
            {
               /* (a) y_b == z_b   =>    o_b == x_b        o_(b+1) == x_b */
               if( b == wordbitstart )
               {
                  CHECK_OKAY( aggregateBinVarsPos(scip, cons, consdata->vars[word], op1vars[b],
                                 naggrvars, infeasible, &redundant) );
               }
               if( b+1 == nextwordbitstart && !(*infeasible) )
               {
                  CHECK_OKAY( aggregateBinVarsPos(scip, cons, consdata->vars[word+1], op1vars[b],
                                 naggrvars, infeasible, &redundant) );
               }
            }
            else
            {
               /* (b) y_b != z_b   =>    o_b == 1 - x_b    o_(b+1) fixed (done in propagation) */
               assert(ovrstate[b+1] == op2state[b]);
               if( b == wordbitstart )
               {
                  CHECK_OKAY( aggregateBinVarsNeg(scip, cons, consdata->vars[word], op1vars[b], 
                                 naggrvars, infeasible, &redundant) );
               }
            }
         }
         else if( ovrfixed && resfixed )
         {
            /* (6) o_b, z_b fixed */
            if( ovrstate[b] == resstate[b] )
            {
               /* (a) o_b == z_b   =>    x_b == y_b        o_(b+1) == y_b */
               CHECK_OKAY( aggregateBinVarsPos(scip, cons, op1vars[b], op2vars[b],
                              naggrvars, infeasible, &redundant) );
               if( b+1 == nextwordbitstart && !(*infeasible) )
               {
                  CHECK_OKAY( aggregateBinVarsPos(scip, cons, consdata->vars[word+1], op2vars[b],
                                 naggrvars, infeasible, &redundant) );
               }
            }
            else
            {
               /* (b) o_b != z_b   =>    x_b == 1 - y_b    o_(b+1) fixed (done in propagation) */
               assert(ovrstate[b+1] == ovrstate[b]);
               CHECK_OKAY( aggregateBinVarsNeg(scip, cons, op1vars[b], op2vars[b], naggrvars, infeasible, &redundant) );
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** propagates domains of variables of an add constraint with the overflow rules:
 *    (4) if o_(b+1) == 0, and z_b == 0, then x_b == y_b == o_b == 0 can be deduced
 *    (5) if o_(b+1) == 0, and one of x_b, y_b, o_b is fixed to 1, the others can be deduced to 0 and z_b to 1
 *    (6) if o_(b+1) == 1, and z_b == 1, then x_b == y_b == o_b == 1 can be deduced
 *    (7) if o_(b+1) == 1, and one of x_b, y_b, o_b is fixed to 0, the others can be deduced to 1 and z_b to 0
 */
static
RETCODE propagateAddOvr(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitarith constraint */
   BITSTATE*        op1state,           /**< current states of operand1 bits */
   BITSTATE*        op2state,           /**< current states of operand2 bits */
   BITSTATE*        ovrstate,           /**< current states of overflow bits */
   BITSTATE*        resstate,           /**< current states of resultant bits */
   VAR**            op1vars,            /**< binary variables representing operand1 */
   VAR**            op2vars,            /**< binary variables representing operand2 */
   VAR**            resvars,            /**< binary variables representing resultant */
   int              b,                  /**< bit number to apply rules on */
   int*             nfixedvars,         /**< pointer to add up the number of fixed variables, or NULL */
   Bool*            infeasible,         /**< pointer to store whether the constraint is infeasible in current bounds */
   Bool*            ovrdeduced          /**< pointer to store whether a new deduction on o_b was made */
   )
{
   CONSDATA* consdata;

   assert(op1state != NULL);
   assert(op2state != NULL);
   assert(ovrstate != NULL);
   assert(resstate != NULL);
   assert(op1vars != NULL);
   assert(op2vars != NULL);
   assert(resvars != NULL);
   assert(infeasible != NULL);
   assert(ovrdeduced != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= b && b < SCIPgetNBitsBitvar(consdata->resultant));

   /* (4) if o_(b+1) == 0, and z_b == 0, then x_b == y_b == o_b == 0 can be deduced */
   /* (5) if o_(b+1) == 0, and one of x_b, y_b, o_b is fixed to 1, the others can be deduced to 0 and z_b to 1 */
   if( ovrstate[b+1] == FIXEDZERO )
   {
      if( resstate[b] == FIXEDZERO )
      {
         if( fixState(op1state, b, FIXEDZERO, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, op1vars[b], 0.0, nfixedvars, infeasible) );
         }
         if( fixState(op2state, b, FIXEDZERO, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, op2vars[b], 0.0, nfixedvars, infeasible) );
         }
         if( fixState(ovrstate, b, FIXEDZERO, infeasible) )
            *ovrdeduced = TRUE;
      }
      else if( op1state[b] == FIXEDONE )
      {
         if( fixState(op2state, b, FIXEDZERO, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, op2vars[b], 0.0, nfixedvars, infeasible) );
         }
         if( fixState(ovrstate, b, FIXEDZERO, infeasible) )
            *ovrdeduced = TRUE;
         if( fixState(resstate, b, FIXEDONE, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, resvars[b], 1.0, nfixedvars, infeasible) );
         }
      }
      else if( op2state[b] == FIXEDONE )
      {
         if( fixState(op1state, b, FIXEDZERO, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, op1vars[b], 0.0, nfixedvars, infeasible) );
         }
         if( fixState(ovrstate, b, FIXEDZERO, infeasible) )
            *ovrdeduced = TRUE;
         if( fixState(resstate, b, FIXEDONE, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, resvars[b], 1.0, nfixedvars, infeasible) );
         }
      }
      else if( ovrstate[b] == FIXEDONE )
      {
         if( fixState(op1state, b, FIXEDZERO, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, op1vars[b], 0.0, nfixedvars, infeasible) );
         }
         if( fixState(op2state, b, FIXEDZERO, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, op2vars[b], 0.0, nfixedvars, infeasible) );
         }
         if( fixState(resstate, b, FIXEDONE, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, resvars[b], 1.0, nfixedvars, infeasible) );
         }
      }
   }

   /* (6) if o_(b+1) == 1, and z_b == 1, then x_b == y_b == o_b == 1 can be deduced */
   /* (7) if o_(b+1) == 1, and one of x_b, y_b, o_b is fixed to 0, the others can be deduced to 1 and z_b to 0 */
   if( ovrstate[b+1] == FIXEDONE )
   {
      if( resstate[b] == FIXEDONE )
      {
         if( fixState(op1state, b, FIXEDONE, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, op1vars[b], 1.0, nfixedvars, infeasible) );
         }
         if( fixState(op2state, b, FIXEDONE, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, op2vars[b], 1.0, nfixedvars, infeasible) );
         }
         if( fixState(ovrstate, b, FIXEDONE, infeasible) )
            *ovrdeduced = TRUE;
      }
      else if( op1state[b] == FIXEDZERO )
      {
         if( fixState(op2state, b, FIXEDONE, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, op2vars[b], 1.0, nfixedvars, infeasible) );
         }
         if( fixState(ovrstate, b, FIXEDONE, infeasible) )
            *ovrdeduced = TRUE;
         if( fixState(resstate, b, FIXEDZERO, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, resvars[b], 0.0, nfixedvars, infeasible) );
         }
      }
      else if( op2state[b] == FIXEDZERO )
      {
         if( fixState(op1state, b, FIXEDONE, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, op1vars[b], 1.0, nfixedvars, infeasible) );
         }
         if( fixState(ovrstate, b, FIXEDONE, infeasible) )
            *ovrdeduced = TRUE; 
         if( fixState(resstate, b, FIXEDZERO, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, resvars[b], 0.0, nfixedvars, infeasible) );
         }
      }
      else if( ovrstate[b] == FIXEDZERO )
      {
         if( fixState(op1state, b, FIXEDONE, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, op1vars[b], 1.0, nfixedvars, infeasible) );
         }
         if( fixState(op2state, b, FIXEDONE, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, op2vars[b], 1.0, nfixedvars, infeasible) );
         }
         if( fixState(resstate, b, FIXEDZERO, infeasible) )
         {
            CHECK_OKAY( deduceBinVar(scip, cons, resvars[b], 0.0, nfixedvars, infeasible) );
         }
      }
   }

   return SCIP_OKAY;
}

/** propagates domains of variables of an add constraint:
 *  A bit addition has the following form:
 *           x_(n-1)  ...  x_0
 *    +      y_(n-1)  ...  y_0
 *    + o_n  o_(n-1)  ...  o_0
 *    ------------------------
 *    =      z_(n-1)  ...  z_0
 *  where each variable is binary. x and y are the operands, o is the overflow and z is the resultant.
 *  The following equalities hold:
 *                o_0 == 0
 *    x_b + y_b + o_b == z_b + 2*o_(b+1)   for all b= 0,...,n-1
 *  The following propagations can be made:
 *    (1) if 3 of x_b, y_b, o_b, and z_b are fixed, the 4th variable and o_(b+1) can be deduced
 *    (2) if x_b + y_b + o_b - z_b <= 1, o_(b+1) == 0 can be deduced
 *    (3) if x_b + y_b + o_b - z_b >= 1, o_(b+1) == 1 can be deduced
 *    (4) if o_(b+1) == 0, and z_b == 0, then x_b == y_b == o_b == 0 can be deduced
 *    (5) if o_(b+1) == 0, and one of x_b, y_b, o_b is fixed to 1, the others can be deduced to 0 and z_b to 1
 *    (6) if o_(b+1) == 1, and z_b == 1, then x_b == y_b == o_b == 1 can be deduced
 *    (7) if o_(b+1) == 1, and one of x_b, y_b, o_b is fixed to 0, the others can be deduced to 1 and z_b to 0
 *  In preprocessing, the following aggregation can be made:
 *    (1) if 2 of x_b, y_b, o_b, and z_b are fixed, one of the other two variables and o_(b+1) can be aggregated
 */
static
RETCODE propagateAdd(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitarith constraint */
   int*             nfixedvars,         /**< pointer to add up the number of fixed variables, or NULL */
   int*             naggrvars,          /**< pointer to add up the number of aggregated variables, or NULL */
   int*             ndelconss,          /**< pointer to add up the number of deleted constraints, or NULL */
   Bool*            infeasible          /**< pointer to store whether the constraint is infeasible in current bounds */
   )
{
   CONSDATA* consdata;
   BITSTATE* op1state;
   BITSTATE* op2state;
   BITSTATE* ovrstate;
   BITSTATE* resstate;
   VAR** op1vars;
   VAR** op2vars;
   VAR** resvars;
   Bool op1fixed;
   Bool op2fixed;
   Bool ovrfixed;
   Bool resfixed;
   Bool ovrdeduced;
   int nop1bits;
   int nop2bits;
   int nresbits;
   int nreswords;
   int bitval;
   int fixsum;
   int w;
   int b;
   int i;

   /* check if the defines for bit states have the expected values */
   assert(FIXEDZERO == -1);
   assert(UNFIXED == 0);
   assert(FIXEDONE == +1);

   assert(infeasible != NULL);
   assert(*infeasible == FALSE);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* get number of bits */
   nop1bits = SCIPgetNBitsBitvar(consdata->operand1);
   nop2bits = SCIPgetNBitsBitvar(consdata->operand2);
   nresbits = SCIPgetNBitsBitvar(consdata->resultant);
   nreswords = SCIPgetNWordsBitvar(consdata->resultant);
   assert(nop1bits <= nresbits);
   assert(nop2bits <= nresbits);

   /* get binary variables for operands and resultant */
   op1vars = SCIPgetBitsBitvar(consdata->operand1);
   op2vars = SCIPgetBitsBitvar(consdata->operand2);
   resvars = SCIPgetBitsBitvar(consdata->resultant);
   assert(op1vars != NULL);
   assert(op2vars != NULL);
   assert(resvars != NULL);

   /* get temporary memory for bit states */
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &op1state, nresbits) );
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &op2state, nresbits) );
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &ovrstate, nresbits+1) );
   CHECK_OKAY( SCIPcaptureBufferArray(scip, &resstate, nresbits) );

   /* initialize bit states for operands and resultant */
   for( b = 0; b < nresbits; ++b )
   {
      if( b >= nop1bits || SCIPvarGetUbLocal(op1vars[b]) < 0.5 )
         op1state[b] = FIXEDZERO;
      else if( SCIPvarGetLbLocal(op1vars[b]) > 0.5 )
         op1state[b] = FIXEDONE;
      else
         op1state[b] = UNFIXED;

      if( b >= nop2bits || SCIPvarGetUbLocal(op2vars[b]) < 0.5 )
         op2state[b] = FIXEDZERO;
      else if( SCIPvarGetLbLocal(op2vars[b]) > 0.5 )
         op2state[b] = FIXEDONE;
      else
         op2state[b] = UNFIXED;

      if( SCIPvarGetUbLocal(resvars[b]) < 0.5 )
         resstate[b] = FIXEDZERO;
      else if( SCIPvarGetLbLocal(resvars[b]) > 0.5 )
         resstate[b] = FIXEDONE;
      else
         resstate[b] = UNFIXED;
   }
   
   /* initialize bit states for overflow bits
    * - only overflows at a word boundary might be fixed
    * - overflow bit 0 is fixed to zero
    */
   clearMemoryArray(ovrstate, nresbits+1); /* set all bits to UNFIXED */
   ovrstate[0] = FIXEDZERO;
   b = 0;
   for( w = 0; w < nreswords; ++w )
   {
      b += SCIPgetWordSizeBitvar(consdata->resultant, w);
      assert(b <= nresbits);
      if( SCIPvarGetUbLocal(consdata->vars[w]) < 0.5 )
         ovrstate[b] = FIXEDZERO;
      else if( SCIPvarGetLbLocal(consdata->vars[w]) > 0.5 )
         ovrstate[b] = FIXEDONE;
   }

   /* apply the propagation rules from least significant bit to most significant bit
    * if an overflow variable was deduced, go back to the lesser significant bits, because propagations
    * (4) to (7) may hold now
    */
   for( b = 0; b < nresbits && !(*infeasible); ++b )
   {
      ovrdeduced = FALSE;

      /*printf("bit %3d: op1=%+d op2=%+d ovr=%+d res=%+d nextovr=%+d\n",
        b, op1state[b], op2state[b], ovrstate[b], resstate[b], ovrstate[b+1]);*/ /*???????????????????*/
      /* (1) if 3 of x_b, y_b, o_b, and z_b are fixed, the 4th variable and o_(b+1) can be deduced */
      fixsum = (op1state[b] == FIXEDONE) + (op2state[b] == FIXEDONE)
         + (ovrstate[b] == FIXEDONE) - (resstate[b] == FIXEDONE);
      op1fixed = (op1state[b] != UNFIXED);
      op2fixed = (op2state[b] != UNFIXED);
      ovrfixed = (ovrstate[b] != UNFIXED);
      resfixed = (resstate[b] != UNFIXED);
      if( op1fixed + op2fixed + ovrfixed + resfixed >= 3 )
      {
         if( !op1fixed )
         {
            assert(op2fixed && ovrfixed && resfixed);

            /* fix x_b to be bit 0 of y_b + o_b - z_b */
            bitval = (fixsum & 0x01);
            if( fixState(op1state, b, bitval*2 - 1, infeasible) )
            {
               CHECK_OKAY( deduceBinVar(scip, cons, op1vars[b], (Real)bitval, nfixedvars, infeasible) );
            }
            fixsum += (op1state[b] == FIXEDONE);
         }
         else if( !op2fixed )
         {
            assert(op1fixed && ovrfixed && resfixed);

            /* fix y_b to be bit 0 of x_b + o_b - z_b */
            bitval = (fixsum & 0x01);
            if( fixState(op2state, b, bitval*2 - 1, infeasible) )
            {
               CHECK_OKAY( deduceBinVar(scip, cons, op2vars[b], (Real)bitval, nfixedvars, infeasible) );
            }
            fixsum += (op2state[b] == FIXEDONE);
         }
         else if( !ovrfixed )
         {
            assert(op1fixed && op2fixed && resfixed);

            /* fix o_b to be bit 0 of x_b + y_b - z_b */
            bitval = (fixsum & 0x01);
            if( fixState(ovrstate, b, bitval*2 - 1, infeasible) )
               ovrdeduced = TRUE;
            fixsum += (ovrstate[b] == FIXEDONE);
         }
         else if( !resfixed )
         {
            assert(op1fixed && op2fixed && ovrfixed);

            /* fix z_b to be bit 0 of x_b + y_b + o_b */
            bitval = (fixsum & 0x01);
            if( fixState(resstate, b, bitval*2 - 1, infeasible) )
            {
               CHECK_OKAY( deduceBinVar(scip, cons, resvars[b], (Real)bitval, nfixedvars, infeasible) );
            }
            fixsum -= (resstate[b] == FIXEDONE);
         }

         /* fix o_(b+1) to be bit 1 of x_b + y_b + o_b - z_b */
         assert((fixsum == 0) || (fixsum == 2));
         (void)fixState(ovrstate, b+1, fixsum - 1, infeasible);
      }
      else
      {
         int sumub;

         /* (2) if x_b + y_b + o_b - z_b <= 1, o_(b+1) == 0 can be deduced */
         sumub = (op1state[b]+2)/2 + (op2state[b]+2)/2 + (ovrstate[b]+2)/2 - (resstate[b]+1)/2;
         if( sumub <= 1 )
            (void)fixState(ovrstate, b+1, FIXEDZERO, infeasible);
         else
         {
            int sumlb;

            /* (3) if x_b + y_b + o_b - z_b >= 1, o_(b+1) == 1 can be deduced */
            sumlb = (op1state[b]+1)/2 + (op2state[b]+1)/2 + (ovrstate[b]+1)/2 - (resstate[b]+2)/2;
            if( sumlb >= 1 )
               (void)fixState(ovrstate, b+1, FIXEDONE, infeasible);
         }

         /* (4) if o_(b+1) == 0, and z_b == 0, then x_b == y_b == o_b == 0 can be deduced */
         /* (5) if o_(b+1) == 0, and one of x_b, y_b, o_b is fixed to 1, the others can be deduced to 0 and z_b to 1 */
         /* (6) if o_(b+1) == 1, and z_b == 1, then x_b == y_b == o_b == 1 can be deduced */
         /* (7) if o_(b+1) == 1, and one of x_b, y_b, o_b is fixed to 0, the others can be deduced to 1 and z_b to 0 */
         CHECK_OKAY( propagateAddOvr(scip, cons, op1state, op2state, ovrstate, resstate, op1vars, op2vars, resvars,
                        b, nfixedvars, infeasible, &ovrdeduced) );
      }
      /*printf("      -> op1=%+d op2=%+d ovr=%+d res=%+d nextovr=%+d   ovrdeduced=%d\n",
        op1state[b], op2state[b], ovrstate[b], resstate[b], ovrstate[b+1], ovrdeduced);*/ /*???????????????????*/

      /* as long as o_b was newly deduced, check the lower bits for further propagations of type (4) to (7) */
      for( i = b-1; ovrdeduced && !(*infeasible); --i )
      {
         assert(i >= 0);
         ovrdeduced = FALSE;
         CHECK_OKAY( propagateAddOvr(scip, cons, op1state, op2state, ovrstate, resstate, op1vars, op2vars, resvars, 
                        i, nfixedvars, infeasible, &ovrdeduced) );
      }
      /*printf("  -----> op1=%+d op2=%+d ovr=%+d res=%+d nextovr=%+d   ovrdeduced=%d\n",
        op1state[b], op2state[b], ovrstate[b], resstate[b], ovrstate[b+1], ovrdeduced);*/ /*???????????????????*/
   }

   /* scan the overflow bits that exist as variables (the word's boundaries) for fixings */
   b = 0;
   for( w = 0; w < nreswords && !(*infeasible); ++w )
   {
      b += SCIPgetWordSizeBitvar(consdata->resultant, w);
      assert(b <= nresbits);
      if( ovrstate[b] == FIXEDZERO && SCIPvarGetUbLocal(consdata->vars[w]) > 0.5 )
      {
         CHECK_OKAY( deduceBinVar(scip, cons, consdata->vars[w], 0.0, nfixedvars, infeasible) );
      }
      else if( ovrstate[b] == FIXEDONE && SCIPvarGetLbLocal(consdata->vars[w]) < 0.5 )
      {
         CHECK_OKAY( deduceBinVar(scip, cons, consdata->vars[w], 1.0, nfixedvars, infeasible) );
      }
   }

   if( !(*infeasible) )
   {
      /* check if we can disable the constraint, because all variables are fixed */
      for( b = 0; b < nop1bits && op1state[b] != UNFIXED; ++b )
      {}
      if( b == nop1bits )
      {
         for( b = 0; b < nop2bits && op2state[b] != UNFIXED; ++b )
         {}
         if( b == nop2bits )
         {
            /* both operands are fixed -> all variables are fixed */
#ifdef DEBUG
            for( b = 0; b < nresbits; ++b )
            {
               assert(ovrstate[b] != UNFIXED);
               assert(resstate[b] != UNFIXED);
            }
            assert(ovrstate[nresbits] != UNFIXED);
#endif
            /* disable the constraint */
            CHECK_OKAY( SCIPdisableConsLocal(scip, cons) );
            if( ndelconss != NULL )
               (*ndelconss)++;
         }
      }
      
      /* in presolving, aggregate variables if possible */
      if( SCIPstage(scip) == SCIP_STAGE_PRESOLVING )
      {
         CHECK_OKAY( aggregateAdd(scip, cons, op1state, op2state, ovrstate, resstate, op1vars, op2vars, resvars,
                        naggrvars, infeasible) );
      }
   }

   /* free temporary memory */
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &resstate) );
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &ovrstate) );
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &op2state) );
   CHECK_OKAY( SCIPreleaseBufferArray(scip, &op1state) );

   return SCIP_OKAY;
}

/** propagates domains of variables of a bitarith constraint */
static
RETCODE propagateCons(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< bitarith constraint */
   int*             nfixedvars,         /**< pointer to add up the number of fixed variables, or NULL */
   int*             naggrvars,          /**< pointer to add up the number of aggregated variables, or NULL */
   int*             nchgbds,            /**< pointer to add up the number of changed bounds, or NULL */
   int*             ndelconss,          /**< pointer to add up the number of deleted constraints, or NULL */
   Bool*            infeasible          /**< pointer to store whether the constraint is infeasible in current bounds */
   )
{
   CONSDATA* consdata;

   assert(infeasible != NULL);

   *infeasible = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check, if the constraint is already propagated */
   if( consdata->propagated )
      return SCIP_OKAY;

   /* mark the constraint propagated */
   consdata->propagated = TRUE;

   /* propagate the bit arithmetic operation */
   switch( consdata->bitarithtype )
   {
   case SCIP_BITARITHTYPE_ADD:
      CHECK_OKAY( propagateAdd(scip, cons, nfixedvars, naggrvars, ndelconss, infeasible) );
      break;

   case SCIP_BITARITHTYPE_SHL:
   case SCIP_BITARITHTYPE_EQ:
   case SCIP_BITARITHTYPE_NOT:
      errorMessage("not implemented yet");
      abort();
      break;

   default:
      errorMessage("invalid bit arithmetic type");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
DECL_CONSFREE(consFreeBitarith)
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
#define consInitBitarith NULL


/** deinitialization method of constraint handler (called when problem solving exits) */
#define consExitBitarith NULL


/** solving start notification method of constraint handler (called when presolving was finished) */
#define consSolstartBitarith NULL


/** frees specific constraint data */
static
DECL_CONSDELETE(consDeleteBitarith)
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
DECL_CONSTRANS(consTransBitarith)
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
   assert(sourcedata->vars == NULL); /* original bitarith constraints must not have any internal variables assigned */
   assert(sourcedata->rows == NULL); /* in original problem, there cannot be LP rows */
   assert(sourcedata->nrows == 0);

   /* create constraint data for target constraint */
   CHECK_OKAY( consdataCreate(scip, &targetdata, sourcedata->bitarithtype, 
                  sourcedata->operand1, sourcedata->operand2, sourcedata->resultant) );

   /* transform the operand and resultant bitvar constraints */
   CHECK_OKAY( consdataTransformOperands(scip, targetdata) );

   /* create target constraint */
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                  SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
                  SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
                  SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), SCIPconsIsRemoveable(sourcecons)) );

   /* create internal variables for constraint */
   assert(SCIPconsIsTransformed(*targetcons));
   CHECK_OKAY( createVars(scip, *targetcons) );

   /* catch variable events */
   CHECK_OKAY( consdataCatchEvents(scip, targetdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
DECL_CONSINITLP(consInitlpBitarith)
{
   int c;

   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* add LP relaxations for all initial constraints */
   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsIsInitial(conss[c]) )
      {
         debugMessage("adding initial bitarith constraint <%s> to LP\n", SCIPconsGetName(conss[c]));
         CHECK_OKAY( addAllCuts(scip, conss[c], 0.0) );
      }
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler */
static
DECL_CONSSEPA(consSepaBitarith)
{
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* step 1: check all useful bitarith constraints for feasibility */
   for( c = 0; c < nusefulconss; ++c )
   {
      debugMessage("separating bitarith constraint <%s>\n", SCIPconsGetName(conss[c]));
      CHECK_OKAY( separateCons(scip, conss[c], result) );
   }

   /* step 2: if no cuts were found and we are in the root node, check remaining bitarith constraints for feasibility */
   if( SCIPgetActDepth(scip) == 0 )
   {
      for( c = nusefulconss; c < nconss && *result == SCIP_DIDNOTFIND; ++c )
      {
         debugMessage("separating bitarith constraint <%s>\n", SCIPconsGetName(conss[c]));
         CHECK_OKAY( separateCons(scip, conss[c], result) );
      }
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpBitarith)
{
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* check for violated constraints
    * LP is processed at current node -> we can add violated bitarith constraints to the LP */

   *result = SCIP_FEASIBLE;

   /* step 1: check all useful bitarith constraints for feasibility */
   for( c = 0; c < nusefulconss; ++c )
   {
      debugMessage("LP enforcing bitarith constraint <%s>\n", SCIPconsGetName(conss[c]));
      CHECK_OKAY( separateCons(scip, conss[c], result) );
   }
   if( *result != SCIP_FEASIBLE )
      return SCIP_OKAY;

   /* step 2: check all obsolete bitarith constraints for feasibility */
   for( c = nusefulconss; c < nconss && *result == SCIP_FEASIBLE; ++c )
   {
      debugMessage("LP enforcing bitarith constraint <%s>\n", SCIPconsGetName(conss[c]));
      CHECK_OKAY( separateCons(scip, conss[c], result) );
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
DECL_CONSENFOPS(consEnfopsBitarith)
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

   /* check all bitarith constraints for feasibility */
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
DECL_CONSCHECK(consCheckBitarith)
{
   Bool violated;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* check all bitarith constraints for feasibility */
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
DECL_CONSPROP(consPropBitarith)
{
   Bool infeasible;
   int nfixedvars;
   int naggrvars;
   int nchgbds;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* propagate all useful bitarith constraints */
   infeasible = FALSE;
   nfixedvars = 0;
   nchgbds = 0;
   for( c = 0; c < nusefulconss && !infeasible; ++c )
   {
      debugMessage("propagating bitarith constraint <%s>\n", SCIPconsGetName(conss[c]));
      CHECK_OKAY( propagateCons(scip, conss[c], &nfixedvars, &naggrvars, &nchgbds, NULL, &infeasible) );
   }

   /* adjust the result */
   if( infeasible )
      *result = SCIP_CUTOFF;
   else if( nfixedvars > 0 || naggrvars > 0 || nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
DECL_CONSPRESOL(consPresolBitarith)
{
   Bool infeasible;
   int nactfixedvars;
   int nactaggrvars;
   int nactchgbds;
   int nactdelconss;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* propagate all bitarith constraints */
   nactfixedvars = 0;
   nactchgbds = 0;
   nactdelconss = 0;
   infeasible = FALSE;
   for( c = 0; c < nconss && !infeasible; ++c )
   {
      CHECK_OKAY( propagateCons(scip, conss[c], &nactfixedvars, &nactaggrvars, &nactchgbds, &nactdelconss, &infeasible) );
   }

   /* adjust the result */
   if( infeasible )
      *result = SCIP_CUTOFF;
   else if( nactfixedvars > 0 || nactaggrvars > 0 || nactchgbds > 0 || nactdelconss > 0 )
   {
      *result = SCIP_SUCCESS;
      (*nfixedvars) += nactfixedvars;
      (*naggrvars) += nactaggrvars;
      (*nchgbds) += nactchgbds;
      (*ndelconss) += nactdelconss;
   }

   return SCIP_OKAY;
}


/** conflict variable resolving method of constraint handler */
#if 0
static
DECL_CONSRESCVAR(consRescvarBitarith)
{
   errorMessage("method of bitarith constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define consRescvarBitarith NULL
#endif


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockBitarith)
{
   CONSDATA* consdata;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* bitarith constraints are equalities: lock all bit variable constraints in both directions */
   CHECK_OKAY( SCIPlockConsVars(scip, consdata->resultant, nlockspos + nlocksneg, nlockspos + nlocksneg) );
   if( consdata->operand1 != NULL )
   {
      CHECK_OKAY( SCIPlockConsVars(scip, consdata->operand1, nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }
   if( consdata->operand2 != NULL )
   {
      CHECK_OKAY( SCIPlockConsVars(scip, consdata->operand2, nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   /* lock all internal variables in both directions */
   for( v = 0; v < consdata->nvars; ++v )
   {
      SCIPvarLock(consdata->vars[v], nlockspos + nlocksneg, nlockspos + nlocksneg);
   }

   return SCIP_OKAY;
}


/** variable rounding unlock method of constraint handler */
static
DECL_CONSUNLOCK(consUnlockBitarith)
{
   CONSDATA* consdata;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* bitarith constraints are equalities: unlock all bit variable constraints in both directions */
   CHECK_OKAY( SCIPunlockConsVars(scip, consdata->resultant, nunlockspos + nunlocksneg, nunlockspos + nunlocksneg) );
   if( consdata->operand1 != NULL )
   {
      CHECK_OKAY( SCIPunlockConsVars(scip, consdata->operand1, nunlockspos + nunlocksneg, nunlockspos + nunlocksneg) );
   }
   if( consdata->operand2 != NULL )
   {
      CHECK_OKAY( SCIPunlockConsVars(scip, consdata->operand2, nunlockspos + nunlocksneg, nunlockspos + nunlocksneg) );
   }

   /* unlock all internal variables in both directions */
   for( v = 0; v < consdata->nvars; ++v )
   {
      SCIPvarUnlock(consdata->vars[v], nunlockspos + nunlocksneg, nunlockspos + nunlocksneg);
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#define consActiveBitarith NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveBitarith NULL


/** constraint enabling notification method of constraint handler */
#define consEnableBitarith NULL


/** constraint disabling notification method of constraint handler */
#define consDisableBitarith NULL




/*
 * bitarith event handler methods
 */

static
DECL_EVENTEXEC(eventExecBitarith)
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

/** creates the handler for bitarith constraints and includes it in SCIP */
RETCODE SCIPincludeConsHdlrBitarith(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create event handler for bound change events */
   CHECK_OKAY( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
                  NULL, NULL, NULL,
                  NULL, eventExecBitarith,
                  NULL) );

   /* create bitarith constraint handler data */
   CHECK_OKAY( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConsHdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_NEEDSCONS,
                  consFreeBitarith, consInitBitarith, consExitBitarith, consSolstartBitarith,
                  consDeleteBitarith, consTransBitarith, consInitlpBitarith,
                  consSepaBitarith, consEnfolpBitarith, consEnfopsBitarith, consCheckBitarith, 
                  consPropBitarith, consPresolBitarith, consRescvarBitarith,
                  consLockBitarith, consUnlockBitarith,
                  consActiveBitarith, consDeactiveBitarith, 
                  consEnableBitarith, consDisableBitarith,
                  conshdlrdata) );

   return SCIP_OKAY;
}

/** creates and captures a bitarith constraint */
RETCODE SCIPcreateConsBitarith(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   BITARITHTYPE     bitarithtype,       /**< type of arithmetic bit operation */
   CONS*            operand1,           /**< bitvar constraint: first (left) operand in operation (x) */
   CONS*            operand2,           /**< bitvar constraint: second (right) operand in operation (y) */
   CONS*            resultant,          /**< bitvar constraint: result of operation (z) */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable          /**< should the constraint be removed from the LP due to aging or cleanup? */
   )
{
   CONSHDLR* conshdlr;
   CONSHDLRDATA* conshdlrdata;
   CONSDATA* consdata;

   /* find the bitarith constraint handler */
   conshdlr = SCIPfindConsHdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("bitarith constraint handler not found");
      return SCIP_PLUGINNOTFOUND;
   }
   
   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   /* check the validity of the operands and resultant with the type of arithmetic operation */
   CHECK_OKAY( validateArith(bitarithtype, operand1, operand2, resultant) );

   /* create constraint data */
   CHECK_OKAY( consdataCreate(scip, &consdata, bitarithtype, operand1, operand2, resultant) );

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
                  local, modifiable, removeable) );

   if( SCIPconsIsTransformed(*cons) )
   {
      /* transform the operand and resultant bitvar constraints */
      CHECK_OKAY( consdataTransformOperands(scip, consdata) );

      /* create internal variables for transformed constraints */
      CHECK_OKAY( createVars(scip, *cons) );

      /* catch variable events */
      CHECK_OKAY( consdataCatchEvents(scip, consdata, conshdlrdata->eventhdlr) );
   }

   return SCIP_OKAY;
}

/** gets number of operands for given arithmetic operation */
int SCIPgetArityBitarith(
   BITARITHTYPE     bitarithtype        /**< type of arithmetic bit operation */
   )
{
   assert(bitarithtype < SCIP_NBITARITHTYPES);

   return bitarithArity[bitarithtype];
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_quadratic.c
 * @ingroup CONSHDLRS
 * @brief  constraint handler for quadratic constraints
 * @author Stefan Vigerske
 * 
 * @todo SCIP might fix variables on +/- infinity; remove them in presolve and take care later
 * @todo constraints that are always feasible w.r.t. local/global bounds should be enabled/disabled (see logicor, setppc)
 * @todo round constraint bounds to integers if all coefficients and variables are (impl.) integer
 * @todo constraints in one variable should be replaced by linear variable or similar
 * @todo recognize and reformulate complementarity constraints (x*y = 0)
 * @todo check if some quadratic terms appear in several constraints and try to simplify (e.g., nous1)
 * @todo skip separation in enfolp if for current LP (check LP id) was already separated
 * @todo don't iterate over hash map, use array additionally
 * @todo watch unbounded variables to enable/disable propagation
 * @todo sort order in bilinvar1/bilinvar2 such that the var which is involved in more terms is in bilinvar1, and use this info propagate and AddLinearReform
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <assert.h>
#include <string.h> /* for strcmp */ 
#include <ctype.h>  /* for isdigit */

#include "scip/cons_quadratic.h"
#include "scip/cons_linear.h"
#include "scip/cons_and.h"
#include "scip/cons_varbound.h"
#include "scip/intervalarith.h"
#include "scip/heur_nlp.h"
#include "scip/nlpi.h"

#ifdef WITH_CONSBRANCHNL
#include "cons_branchnonlinear.h"
#endif

/* #define USE_RELAXBRANCH */

/* constraint handler properties */
#define CONSHDLR_NAME          "quadratic"
#define CONSHDLR_DESC          "quadratic constraints of the form lhs <= b' x + x' A x <= rhs"
#define CONSHDLR_SEPAPRIORITY        10 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       -50 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -4000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             2 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            10 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */


/*
 * Data structures
 */

/** Eventdata for variable bound change events. */
struct SCIP_EventData
{
   SCIP_CONSDATA*        consdata;           /**< the constraint data */
   int                   varidx;             /**< the index of the variable which bound change is catched, positive for linear variables, negative for quadratic variables */
};

/** Data of a quadratic constraint. */
struct SCIP_ConsData
{
   SCIP_Real             lhs;                /**< left hand side of constraint */
   SCIP_Real             rhs;                /**< right hand side of constraint */

   int                   nlinvars;           /**< number of linear variables */
   SCIP_VAR**            linvars;            /**< linear variables */
   SCIP_Real*            lincoefs;           /**< coefficients of linear variables */

   int                   nquadvars;          /**< number of variables in quadratic terms */
   SCIP_VAR**            quadvars;           /**< array with variables in quadratic terms */
   SCIP_Real*            quadlincoefs;       /**< array with linear coefficients of quadratic variables */
   SCIP_Real*            quadsqrcoefs;       /**< array with coefficients of square terms of quadratic variables */
   int*                  nadjbilin;          /**< number of bilinear terms where the variable is involved */
   int**                 adjbilin;           /**< array with indices of bilinear terms in which variable is involved */

   int                   nbilinterms;        /**< number of bilinear terms */
   SCIP_VAR**            bilinvars1;         /**< array with first variables in bilinear term */
   SCIP_VAR**            bilinvars2;         /**< array with second variables in bilinear term */
   SCIP_Real*            bilincoefs;         /**< array with coefficients of bilinear term */

   unsigned int          isconvex:1;         /**< is quadratic function is convex ? */
   unsigned int          isconcave:1;        /**< is quadratic function is concave ? */
   unsigned int          isremovedfixings:1; /**< did we removed fixed/aggr/multiaggr variables ? */
   unsigned int          ispropagated:1;     /**< was the constraint propagated with respect to the current bounds ? */
   unsigned int          ispresolved:1;      /**< did we checked for possibilities of upgrading or implicit integer variables ? */

   SCIP_Real             lhsviol;            /**< violation of lower bound by current solution (used temporarily inside constraint handler) */
   SCIP_Real             rhsviol;            /**< violation of lower bound by current solution (used temporarily inside constraint handler) */

   SCIP_EVENTDATA*       linbndchgeventdata; /**< eventdata for bound change of linear variable */
   SCIP_EVENTDATA*       quadbndchgeventdata;/**< eventdata for bound change on quadratic variable */
   SCIP_INTERVAL*        linrange;           /**< range of each linear term */
   SCIP_INTERVAL         quadrange;          /**< range of quadratic term as used in isIntervalFeasible */
   SCIP_INTERVAL*        quadrangevar;       /**< range of quadratic term except one variable as used in propagation */
   SCIP_INTERVAL*        bilinrange;         /**< range of bilinear term as used in propagation */
};

/** quadratic term as used during presolve */
struct PresolveQuadTerm
{
   SCIP_Real             lincoefs;           /**< linear coefficient of a variable */
   SCIP_Real             sqrcoefs;           /**< square coefficient of a variable */
   SCIP_HASHMAP*         bilin;              /**< bilinear terms involving a variable: mapping from SCIP_VAR* to PRESOLVEBILINITEM* */
   int                   component;          /**< component(=block) number in block separable form */
};

typedef struct PresolveQuadTerm PRESOLVEQUADTERM; /* quadratic term used during presolving */

/** bilinear term as used inside PresolveQuadTerm */
struct PresolveBilinItem
{
   SCIP_Real             coef;               /**< coefficient of bilinear term */
   int                   bilinidx;           /**< index in bilin array */
};
typedef struct PresolveBilinItem PRESOLVEBILINITEM; /* bilinear term used during presolving */

#ifndef WITH_CONSBRANCHNL
/** Stores information about the infeasibility assigned to a variable */
typedef struct VarInfeasibility VARINFEASIBILITY; /**< Infeasibility of a variable. */
struct VarInfeasibility
{
   SCIP_Real             min;            /**< minimal infeasibility */
   SCIP_Real             max;            /**< maximal infeasibility */
   SCIP_Real             sum;            /**< sum of infeasibilities */
   VARINFEASIBILITY*     next;           /**< pointer to next vars infeasibility */
};
#endif

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             replacesqrbinary;          /**< were squares of binary variables replaced by the variable itself ? */
   int                   replacebinaryprodlength;   /**< length of linear term which when multiplied with a binary variable is replaced by an auxiliary variable and an equivalent linear formulation */
   int                   empathy4and;               /**< how much empathy we have for using the AND constraint handler: 0 avoid always; 1 use sometimes; 2 use as often as possible */
   SCIP_Bool             disaggregation;            /**< should we disaggregate block separable quadratic constraints ? */
   SCIP_Real             mincutefficacy;            /**< minimal efficacy of a cut in order to add it to relaxation */
   SCIP_Bool             doscaling;                 /**< should constraints be scaled in the feasibility check ? */
   SCIP_Bool             fastpropagate;             /**< should a faster but maybe less effective propagation be used ? */
   SCIP_Real             defaultbound;              /**< a bound to set for variables that are unbounded and in a nonconvex term after presolve */
   SCIP_Real             cutmaxrange;               /**< maximal range (maximal coef / minimal coef) of a cut in order to be added to LP */
   SCIP_Bool             linearizenlpsol;           /**< whether convex quadratic constraints should be linearized in a solution found by the NLP or RENSNL heuristic */

   SCIP_HEUR*            nlpheur;                   /**< a pointer to the NLP heuristic, if available */
   SCIP_HEUR*            rensnlheur;                /**< a pointer to the RENSNL heuristic, if available */
   SCIP_EVENTHDLR*       eventhdlr;                 /**< our handler for variable bound change events */
   int                   newsoleventfilterpos;      /**< filter position of new solution event handler, if catched */

#ifdef WITH_CONSBRANCHNL
   SCIP_CONSHDLR*        branchnl;                  /**< constraint handler for branching on nonlinear variables */
#else
#ifndef USE_RELAXBRANCH
   SCIP_HASHMAP*         branchcand;                /**< branching candidates */
   VARINFEASIBILITY*     varinfeas;                 /**< list of variable infeasibilities */
   
   char                  strategy;                  /**< branching strategy */
#endif
   SCIP_Real             mindistbrpointtobound;     /**< minimal (fractional) distance of branching point to bound */
#endif
   
   SCIP_QUADCONSUPGRADE** quadconsupgrades;         /**< quadratic constraint upgrade methods for specializing quadratic constraints */
   int                   quadconsupgradessize;      /**< size of quadconsupgrade array */
   int                   nquadconsupgrades;         /**< number of quadratic constraint upgrade methods */
#ifdef USECLOCK   
   SCIP_CLOCK*           clock1;
   SCIP_CLOCK*           clock2;
   SCIP_CLOCK*           clock3;
#endif
};

/** quadratic constraint update method */
struct SCIP_QuadConsUpgrade
{
   SCIP_DECL_QUADCONSUPGD((*quadconsupgd));    /**< method to call for upgrading quadratic constraint */
   int                     priority;           /**< priority of upgrading method */
   SCIP_Bool               active;             /**< is upgrading enabled */
};

/** data structure for tokenizing strings */
struct SCIP_Tokenizer
{
   const char* strbuf;                       /**< string to be tokenized */
   char* token;                              /**< current token */
   char* tokenbuf;                           /**< token buffer */
   int strpos;                               /**< current position in string */
};
typedef struct SCIP_Tokenizer SCIP_TOKENIZER;     /**< data structure for tokenizing a string */

/** enum type for constraint sense */
enum CipSense
{
   CIP_SENSE_NOTHING, CIP_SENSE_LE, CIP_SENSE_GE, CIP_SENSE_EQ
};
typedef enum CipSense CIPSENSE;                   /**< enum type for constraint sense */

/*
 * Local methods
 */

/* local methods for parsing quadratic constraints */
static const char delimchars[] = " \f\n\r\t\v";
static const char tokenchars[] = "-+<>=^";

/** returns whether the given character is a token delimiter */
static
SCIP_Bool isDelimChar(
   char                  c                   /**< input character */
   )
{
   return (c == '\0') || (strchr(delimchars, c) != NULL);
}

/** returns whether the given character is a single token */
static
SCIP_Bool isTokenChar(
   char                  c                   /**< input character */
   )
{
   return (strchr(tokenchars, c) != NULL);
}

/** returns whether the current token is an equation sense */
static
SCIP_Bool isSense(
   SCIP_TOKENIZER*       tokenizer,          /**< LP reading data */
   CIPSENSE*             sense               /**< pointer to store the equation sense */
   )
{
   assert(tokenizer != NULL);
   assert(sense != NULL);

   if( strcmp(tokenizer->token, "<=") == 0 )
   {
      *sense = CIP_SENSE_LE;
      return TRUE;
   }
   else if( strcmp(tokenizer->token, ">=") == 0 )
   {
      *sense = CIP_SENSE_GE;
      return TRUE;
   }
   else if( strcmp(tokenizer->token, "==") == 0 )
   {
      *sense = CIP_SENSE_EQ;
      return TRUE;
   }

   return FALSE;
}

/** returns whether the current token is a sign */
static
SCIP_Bool isSign(
   SCIP_TOKENIZER*       tokenizer,          /**< LP reading data */
   int*                  sign                /**< pointer to update the sign */
   )
{
   assert(tokenizer != NULL);
   assert(sign != NULL);
   assert(*sign == +1 || *sign == -1);

   if( tokenizer->token[1] == '\0' )
   {
      if( *tokenizer->token == '+' )
         return TRUE;
      else if( *tokenizer->token == '-' )
      {
         *sign *= -1;
         return TRUE;
      }
   }

   return FALSE;
}

/** returns whether the current token is a value */
static
SCIP_Bool isValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TOKENIZER*       lpinput,            /**< LP reading data */
   SCIP_Real*            value               /**< pointer to store the value (unchanged, if token is no value) */
   )
{
   assert(lpinput != NULL);
   assert(value != NULL);

   if( strcasecmp(lpinput->token, "INFINITY") == 0 || strcasecmp(lpinput->token, "INF") == 0 )
   {
      *value = SCIPinfinity(scip);
      return TRUE;
   }
   else
   {
      double val;
      char* endptr;

      val = strtod(lpinput->token, &endptr);
      if( endptr != lpinput->token && *endptr == '\0' )
      {
         *value = val;
         return TRUE;
      }
   }

   return FALSE;
}

/** returns whether the current character is member of a value string */
static
SCIP_Bool isValueChar(
   char                  c,                  /**< input character */
   char                  nextc,              /**< next input character */
   SCIP_Bool*            hasdot              /**< pointer to update the dot flag */
   )
{
   assert(hasdot != NULL);

   if( isdigit(c) )
      return TRUE;
   else if( !(*hasdot) && (c == '.') && isdigit(nextc) )
   {
      *hasdot = TRUE;
      return TRUE;
   }

   return FALSE;
}

/** reads the next token from the string into the token buffer; returns whether a token was read */
static
SCIP_Bool getNextToken(
   SCIP_TOKENIZER*       tokenizer           /**< string tokenizer */
   )
{
   SCIP_Bool hasdot;
   int tokenlen;
   const char* buf;
   
   /* skip delimiters */
   buf = tokenizer->strbuf;
   while( isDelimChar(buf[tokenizer->strpos]) )
   {
      if( buf[tokenizer->strpos] == '\0' )
      {
         SCIPdebugMessage("end of string\n");
         return FALSE;
      }
      
      tokenizer->strpos++;
   }

   assert(!isDelimChar(buf[tokenizer->strpos]));

   /* check if the token is a value */
   hasdot = FALSE;
   if( isValueChar(buf[tokenizer->strpos], buf[tokenizer->strpos+1], &hasdot) )
   {
      /* read value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < SCIP_MAXSTRLEN);
         assert(!isDelimChar(buf[tokenizer->strpos]));
         tokenizer->token[tokenlen] = buf[tokenizer->strpos];
         tokenlen++;
         tokenizer->strpos++;
      }
      while( isValueChar(buf[tokenizer->strpos], buf[tokenizer->strpos+1], &hasdot) );
   }
   else
   {
      /* read non-value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < SCIP_MAXSTRLEN);
         tokenizer->token[tokenlen] = buf[tokenizer->strpos];
         tokenlen++;
         tokenizer->strpos++;
         if( tokenlen == 1 && isTokenChar(tokenizer->token[0]) )
            break;
      }
      while( !isDelimChar(buf[tokenizer->strpos]) && !isTokenChar(buf[tokenizer->strpos]) );

      /* if the token is an equation sense '<', '>', or '=', skip a following '='
       * if the token is an equality token '=' and the next character is a '<' or '>', replace the token by the inequality sense
       * if the token is an ^ and the next character is a '2', read the next character 
       */
      if( tokenlen == 1
         && (tokenizer->token[0] == '<' || tokenizer->token[0] == '>' || tokenizer->token[0] == '=')
         && buf[tokenizer->strpos] == '=' )
      {
         tokenizer->token[tokenlen] = buf[tokenizer->strpos];
         tokenlen++;
         tokenizer->strpos++;
      }
      else if( tokenlen == 1 && tokenizer->token[0] == '=' && (buf[tokenizer->strpos] == '<' || buf[tokenizer->strpos] == '>') )
      {
         tokenlen++;
         tokenizer->token[1] = tokenizer->token[0];
         tokenizer->token[0] = buf[tokenizer->strpos];
         tokenizer->token[1] = buf[tokenizer->strpos];
         tokenizer->strpos++;
      }
      else if( tokenlen == 1 && tokenizer->token[0] == '^' && buf[tokenizer->strpos] == '2' )
      {
         tokenizer->token[1] = buf[tokenizer->strpos];
         tokenlen++;
         tokenizer->strpos++;
      }
   }
   assert(tokenlen < SCIP_MAXSTRLEN);
   tokenizer->token[tokenlen] = '\0';

   SCIPdebugMessage("read token: '%s'\n", tokenizer->token);

   return TRUE;
}

#ifndef WITH_CONSBRANCHNL
#ifndef USE_RELAXBRANCH
/** clears list of branching candidates */
static
SCIP_RETCODE clearBranchingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler  */
   )
{
   SCIP_CONSHDLRDATA* data;
   VARINFEASIBILITY* v;
   VARINFEASIBILITY* w;
  
   assert(scip != NULL);
   assert(conshdlr != NULL);

   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);
   assert(data->branchcand != NULL);
   
   v = data->varinfeas;
   while( v != NULL )
   {
      w = v->next;
      SCIPfreeBlockMemory(scip, &v);
      v = w;
   }
   data->varinfeas = NULL;
   
   SCIP_CALL( SCIPhashmapRemoveAll(data->branchcand) );
   
   return SCIP_OKAY;
}

/** determines branching point for a variable */
static
SCIP_RETCODE selectBranchingPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for branching on nonlinear variables */
   SCIP_VAR*             var,                /**< branching variables */
   SCIP_Real*            leftub,             /**< buffer to store new upper bound of variable in left branch */
   SCIP_Real*            rightlb             /**< buffer to store new lower bound of variable in right branch */
   )
{
   SCIP_CONSHDLRDATA* data;
   SCIP_Real branchpoint;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(var != NULL);
   assert(leftub != NULL);
   assert(rightlb != NULL);
     
   assert(scip != NULL);
   assert(conshdlr != NULL);
  
   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);

   branchpoint = SCIPgetVarSol(scip, var);
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   if( SCIPisInfinity(scip, branchpoint) )
   {
      if( SCIPisPositive(scip, lb) )
         branchpoint = lb + 1000;
      else
         branchpoint = 0.0;
   }
   else if( SCIPisInfinity(scip, -branchpoint) )
   {
      if( SCIPisNegative(scip, ub) )
         branchpoint = ub - 1000;
      else
         branchpoint = 0.0;
   }

   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
   {
      if( !SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub) )
      {
         /* branch on value of LP solution
          * if it is too close to the bounds, move more into the middle of the interval */
         if( branchpoint < (1-data->mindistbrpointtobound) * lb + data->mindistbrpointtobound * ub )
            branchpoint = (1-data->mindistbrpointtobound) * lb + data->mindistbrpointtobound * ub;
         else if( branchpoint > data->mindistbrpointtobound * lb + (1-data->mindistbrpointtobound) * ub )
            branchpoint = data->mindistbrpointtobound * lb + (1-data->mindistbrpointtobound) * ub;

         /* for very tiny intervals we set it into the middle */
         if( !SCIPisGT(scip, branchpoint, lb) || !SCIPisLT(scip, branchpoint, ub) )
            branchpoint = (lb+ub) * .5;
      }
      else if( !SCIPisLT(scip, lb, branchpoint) )
      {
         assert(SCIPisInfinity(scip, ub));
         branchpoint = lb + MAX(0.5*ABS(lb), 1000);
      }
      else if( !SCIPisGT(scip, ub, branchpoint) )
      {
         assert(SCIPisInfinity(scip, -lb));
         branchpoint = ub - MAX(0.5*ABS(ub), 1000);
      }

      *leftub = *rightlb = branchpoint;
   }
   else
   {
      if( branchpoint > ub )
         branchpoint = ub;
      else if( branchpoint < lb )
         branchpoint = lb;
      if( SCIPisIntegral(scip, branchpoint) )
      {
         if( branchpoint < .5*(lb+ub) )
            branchpoint += .5;
         else
            branchpoint -= .5;
      }
      *rightlb = SCIPceil(scip, branchpoint);
      *leftub  = SCIPfloor(scip, branchpoint);
   }
   
   return SCIP_OKAY;
}

/** selects a branching variable and decides about branching point */
static
SCIP_RETCODE selectBranchingVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for branching on nonlinear variables */
   SCIP_VAR**            var,                /**< buffer to store branching variable */
   SCIP_Real*            leftub,             /**< buffer to store new upper bound of branching variable in  left branch */
   SCIP_Real*            rightlb             /**< buffer to store new lower bound of branching varialbe in right branch */
   )
{  /*lint --e{777} */
   SCIP_CONSHDLRDATA*   data;
   int                  listidx;
   SCIP_HASHMAPLIST*    candlist;
   SCIP_VAR*            cand;
   VARINFEASIBILITY*    infeas;
   SCIP_Real            candleftub, candrightlb;
   SCIP_Real            deltaminus, deltaplus;
   SCIP_Real            pscostdown, pscostup;
   SCIP_Real            score, bestscore;
  
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(var != NULL);
   assert(leftub != NULL);
   assert(rightlb != NULL);
  
   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);
   assert(data->branchcand != NULL);
   
   *var = NULL;
   bestscore = -1.0;
   
   for( listidx = 0; listidx < SCIPhashmapGetNLists(data->branchcand); ++listidx )
   {
      for( candlist = SCIPhashmapGetList(data->branchcand, listidx); candlist != NULL; candlist = SCIPhashmapListGetNext(candlist) )
      {
         cand = (SCIP_VAR*) SCIPhashmapListGetOrigin(candlist);
         infeas = (VARINFEASIBILITY*) SCIPhashmapListGetImage(candlist);
         
         SCIP_CALL( selectBranchingPoint(scip, conshdlr, cand, &candleftub, &candrightlb) );
         assert(candleftub == candrightlb || SCIPvarGetType(cand) <= SCIP_VARTYPE_IMPLINT);

         switch( data->strategy )
         {
            case 'b':
               if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)) )
                  deltaminus = SCIPisInfinity(scip, infeas->max) ? SCIPinfinity(scip) : 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max;
               else
                  deltaminus = candleftub - SCIPvarGetLbLocal(cand);
              
               if( SCIPisInfinity(scip,  SCIPvarGetUbLocal(cand)) )
                  deltaplus  = SCIPisInfinity(scip, infeas->max) ? SCIPinfinity(scip) : 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max;
               else
                  deltaplus  = SCIPvarGetUbLocal(cand) - candrightlb;
              
               break;
              
            case 'r':
               if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)) )
                  deltaplus  = SCIPisInfinity(scip, infeas->max) ? SCIPinfinity(scip) : 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max;
               else
                  deltaplus  = candleftub - SCIPvarGetLbLocal(cand);
            
               if( SCIPisInfinity(scip,  SCIPvarGetUbLocal(cand)) )
                  deltaminus = SCIPisInfinity(scip, infeas->max) ? SCIPinfinity(scip) : 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max;
               else
                  deltaminus = SCIPvarGetUbLocal(cand) - candrightlb;
            
               break;

            case 'i':
               deltaminus = deltaplus = 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max;
               break;

            default :
               SCIPerrorMessage("branching strategy %c unknown\n", data->strategy);
               return SCIP_ERROR;
         }

         if( SCIPisInfinity(scip, deltaminus) || SCIPisInfinity(scip, deltaplus) )
            score = SCIPinfinity(scip);
         else
         {
            pscostdown = SCIPgetVarPseudocost(scip, cand, -deltaminus);
            pscostup = SCIPgetVarPseudocost(scip, cand,  deltaplus);
            score = SCIPgetBranchScore(scip, cand, pscostdown, pscostup);
         }
         SCIPdebugMessage("branching score variable %s = %g; \tinfeas = %g; \ttype=%d  bestscore=%g\n", 
            SCIPvarGetName(cand), score, 0.1 * infeas->sum + 0.8 * infeas->min + 1.3 * infeas->max, SCIPvarGetType(cand), bestscore);

         if( SCIPisSumGT(scip, score, bestscore) )
         {
            bestscore = score;
            *var = cand;
            *leftub = candleftub;
            *rightlb = candrightlb;
         }
         else if( SCIPisSumEQ(scip, score, bestscore) && !(SCIPisInfinity(scip, -SCIPvarGetLbLocal(*var)) && SCIPisInfinity(scip, SCIPvarGetUbLocal(*var))) )
         { /* if best candidate so far is bounded or unbounded at atmost one side, maybe take new candidate */
            if( (SCIPisInfinity(scip, -SCIPvarGetLbLocal(cand)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(cand))) &&
                (SCIPisInfinity(scip, -SCIPvarGetLbLocal(*var)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(*var))) )
            { /* if both variables are unbounded but one of them is bounded on one side, take the one with the larger bound on this side 
               * (hope that this avoids branching on always the same variable) */
               if( SCIPvarGetUbLocal(cand) > SCIPvarGetUbLocal(*var) ||
                   SCIPvarGetLbLocal(cand) < SCIPvarGetLbLocal(*var) )
               {
                  *var = cand;
                  *leftub = candleftub;
                  *rightlb = candrightlb;
               }
            }
            else if( SCIPvarGetType(*var) == SCIPvarGetType(cand) )
            { /* if both have the same type, take the one with larger diameter */
               if( SCIPisLT(scip, SCIPvarGetUbLocal(*var) - SCIPvarGetLbLocal(*var), SCIPvarGetUbLocal(cand) - SCIPvarGetLbLocal(cand)) )
               {
                  *var = cand;
                  *leftub = candleftub;
                  *rightlb = candrightlb;
               }
            }
            else if( SCIPvarGetType(*var) > SCIPvarGetType(cand) )
            { /* take the one with better type ("more discrete") */
               *var = cand;
               *leftub = candleftub;
               *rightlb = candrightlb;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** finds a branching variable and does branching */
static
SCIP_RETCODE enforceByBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_RESULT*          result              /**< buffer where to store result */
	)
{  /*lint --e{777} */
   SCIP_CONSHDLRDATA*   data;
   SCIP_VAR*            brvar = NULL;
   SCIP_Real            leftub; 
   SCIP_Real            rightlb;
   SCIP_Real            leftobjest;
   SCIP_Real            rightobjest;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);

   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);
   
   if( data->varinfeas == NULL )
   { /* have no candidates for branching */
      return SCIP_OKAY;
   }

   leftub = 0.0;
   rightlb = 0.0;

   SCIP_CALL( selectBranchingVariable(scip, conshdlr, &brvar, &leftub, &rightlb) );
   SCIP_CALL( clearBranchingCandidates(scip, conshdlr) );
   
   if( brvar == NULL )
   {
      SCIPwarningMessage("branching variable selection failed to select a variable\n");
      return SCIP_OKAY;
   }
   
   leftobjest = SCIPcalcChildEstimate(scip, brvar, leftub);
   rightobjest = (leftub != rightlb) ? SCIPcalcChildEstimate(scip, brvar, rightlb) : leftobjest;
   if( leftobjest  > SCIPinfinity(scip) )
      leftobjest = SCIPinfinity(scip)/5.;
   if( rightobjest > SCIPinfinity(scip) )
      rightobjest = leftobjest;

   if( SCIPvarGetStatus(brvar) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_NODE* node;
      SCIP_CONS* cons;
      SCIP_Real val;
      val = 1.0;
      SCIPdebugMessage("branching on multiaggregated variable %s: new intervals: [%g, %g] [%g, %g]\n", 
         SCIPvarGetName(brvar), SCIPvarGetLbLocal(brvar), leftub, rightlb, SCIPvarGetUbLocal(brvar));

      SCIP_CALL( SCIPcreateChild(scip, &node, 0.0, leftobjest) );
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "branch", 1, &brvar, &val, SCIPvarGetLbLocal(brvar), 
            leftub, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      SCIP_CALL( SCIPcreateChild(scip, &node, 0.0, rightobjest) );
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "branch", 1, &brvar, &val, rightlb, SCIPvarGetUbLocal(brvar), 
            TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }
   else
   {
      if( SCIPvarGetType(brvar) != SCIP_VARTYPE_CONTINUOUS )
      {
         SCIPdebugMessage("branching on discrete variable %s\n", SCIPvarGetName(brvar));
         SCIP_CALL( SCIPbranchVar(scip, brvar, NULL, NULL, NULL) );
      }
      else
      {
         SCIP_NODE* node;
         SCIPdebugMessage("branching on continuous variable %s: new intervals: [%g, %g] [%g, %g]\n", 
            SCIPvarGetName(brvar), SCIPvarGetLbLocal(brvar), leftub, rightlb, SCIPvarGetUbLocal(brvar));

         SCIP_CALL( SCIPcreateChild(scip, &node, 0.0, leftobjest) );
         SCIP_CALL( SCIPchgVarUbNode(scip, node, brvar, leftub) );

         SCIP_CALL( SCIPcreateChild(scip, &node, 0.0, rightobjest) );
         SCIP_CALL( SCIPchgVarLbNode(scip, node, brvar, rightlb) );
      }
   }

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** updates or initializes the infeasibility of a variable
 * 
 *  If called the first time for some variable, then this variable is added to the list of branching candidates.
 */
static
SCIP_RETCODE updateVarInfeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             varinfeasibility    /**< infeasibility of variable */
   )
{
   SCIP_CONSHDLRDATA*  data;
   VARINFEASIBILITY*   varinfeas;
  
   SCIPdebugMessage("register infeasibility %g for variable %s  [%g, %g]\n", 
      varinfeasibility, SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   /* assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0); */
   assert(varinfeasibility >= 0.0);
   assert(!SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));
  
   data = SCIPconshdlrGetData(conshdlr);
   assert(data != NULL);
   assert(data->branchcand != NULL);
   
   varinfeas = (VARINFEASIBILITY*)SCIPhashmapGetImage(data->branchcand, (void*)var);
   
   if( varinfeas == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &varinfeas) );
      varinfeas->min = varinfeasibility;
      varinfeas->max = varinfeasibility;
      varinfeas->sum = varinfeasibility;
      varinfeas->next = data->varinfeas;
      data->varinfeas = varinfeas;
      
      SCIP_CALL( SCIPhashmapInsert(data->branchcand, (void*)var, varinfeas) );
   }
   else
   {
      varinfeas->sum += varinfeasibility;
      if( varinfeasibility < varinfeas->min )
         varinfeas->min = varinfeasibility;
      if( varinfeasibility > varinfeas->max )
         varinfeas->max = varinfeasibility;
   }
   
   return SCIP_OKAY;
}
#endif
#endif

/** translate from one value of infinity to another
 * 
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) (val >= infty1 ? infty2 : val)

#ifndef WITH_LAPACK

/** if WITH_LAPACK not set, but WITH_IPOPT, then we can use Lapack from Ipopt and also get the Fortran naming convention from it */
#ifdef WITH_IPOPT
#define WITH_LAPACK
#include "IpoptConfig.h"
#endif

#else

/** if WITH_LAPACK is set, then also F77_FUNC should be set, otherwise we try a default that works on common systems */
#ifndef F77_FUNC
#warning "do not know about fortran naming convention for using Lapack; please consider defining F77_FUNC"
/* this is compiler and machine dependent; the following just assumes a Linux/gcc system */
#define F77_FUNC(name,NAME) name ## _
/* #define F77_FUNC_(name,NAME) name ## _ */
#endif

#endif /* ifndef/else WITH_LAPACK */

#ifdef WITH_LAPACK

/** LAPACK Fortran subroutine DSYEV */
void F77_FUNC(dsyev,DSYEV)(
   char*                 jobz,               /**< 'N' to compute eigenvalues only, 'V' to compute eigenvalues and eigenvectors */
   char*                 uplo,               /**< 'U' if upper triangle of A is stored, 'L' if lower triangle of A is stored */
   int*                  n,                  /**< dimension */
   double*               A,                  /**< matrix A on entry; orthonormal eigenvectors on exit, if jobz == 'V' and info == 0; if jobz == 'N', then the matrix data is destroyed */
   int*                  ldA,                /**< leading dimension, probably equal to n */ 
   double*               W,                  /**< buffer for the eigenvalues in ascending order */
   double*               WORK,               /**< workspace array */
   int*                  LWORK,              /**< length of WORK; if LWORK = -1, then the optimal workspace size is calculated and returned in WORK(1) */
   int*                  info                /**< == 0: successful exit; < 0: illegal argument at given position; > 0: failed to converge */
);

static
SCIP_RETCODE LapackDsyev(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             computeeigenvectors,/**< should also eigenvectors should be computed ? */
   int                   N,                  /**< dimension */
   SCIP_Real*            a,                  /**< matrix data on input (size N*N); eigenvectors on output if computeeigenvectors == TRUE */
   SCIP_Real*            w                   /**< buffer to store eigenvalues (size N) */
   )
{
   int     INFO;
   char    JOBZ = computeeigenvectors ? 'V' : 'N';
   char    UPLO = 'L';
   int     LDA  = N;
   double* WORK = NULL;
   int     LWORK;
   double  WORK_PROBE;
   int     i;

   /* First we find out how large LWORK should be */
   LWORK = -1;
   F77_FUNC(dsyev,DSYEV)(&JOBZ, &UPLO, &N, a, &LDA, w, &WORK_PROBE, &LWORK, &INFO);
   if( INFO != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEV. INFO = %d\n", INFO);
      return SCIP_ERROR;
   }

   LWORK = (int) WORK_PROBE;
   assert(LWORK > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &WORK, LWORK) );
   for( i = 0; i < LWORK; ++i )
      WORK[i] = i;
   F77_FUNC(dsyev,DSYEV)(&JOBZ, &UPLO, &N, a, &LDA, w, WORK, &LWORK, &INFO);
   SCIPfreeBufferArray(scip, &WORK);
   if( INFO != 0 )
   {
       SCIPerrorMessage("There was an error when calling DSYEV. INFO = %d\n", INFO);
       return SCIP_ERROR;
   }

   return SCIP_OKAY;
}
#endif

static
SCIP_DECL_EVENTEXEC(processVarEvent)
{
   SCIP_CONSDATA* consdata;
   
   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   consdata = eventdata->consdata;
   assert(consdata != NULL);

   if( SCIPeventGetType(event) & SCIP_EVENTTYPE_VARFIXED )
   {
      consdata->isremovedfixings = FALSE;
      /* return SCIP_OKAY; */
   }

   if( consdata->ispropagated )
   { /* if we have a bound tightening, we might not have propagated bounds anymore */
      if( SCIPeventGetType(event) == SCIP_EVENTTYPE_UBTIGHTENED )
      {
         consdata->ispropagated =
            eventdata->varidx >= 0 &&
            (SCIPisInfinity(scip, -consdata->lhs) || consdata->lincoefs[eventdata->varidx] < 0) &&
            (SCIPisInfinity(scip,  consdata->rhs) || consdata->lincoefs[eventdata->varidx] > 0);
      }
      else if( SCIPeventGetType(event) == SCIP_EVENTTYPE_LBTIGHTENED )
      {
         consdata->ispropagated =
            eventdata->varidx >= 0 &&
            (SCIPisInfinity(scip, -consdata->lhs) || consdata->lincoefs[eventdata->varidx] > 0) &&
            (SCIPisInfinity(scip,  consdata->rhs) || consdata->lincoefs[eventdata->varidx] < 0);
      }
   }

   if( eventdata->varidx >= 0 )
   { /* make linrange[varidx] invalid */
      SCIPintervalSetEmpty(&consdata->linrange[eventdata->varidx]);
   }
   else
   { /* make quadrange, quadrangevar[-varidx-1], and some bilinrange invalid */
      assert(consdata->nquadvars > 0);
      assert(-eventdata->varidx-1 < consdata->nquadvars);
      SCIPintervalSetEmpty(&consdata->quadrange);
      if( !SCIPintervalIsEmpty(consdata->quadrangevar[-eventdata->varidx-1]) )
      {
         int i;
         SCIPintervalSetEmpty(&consdata->quadrangevar[-eventdata->varidx-1]);
         for( i = 0; i < consdata->nadjbilin[-eventdata->varidx-1]; ++i )
            SCIPintervalSetEmpty(&consdata->bilinrange[consdata->adjbilin[-eventdata->varidx-1][i]]);
      }
   }

   return SCIP_OKAY;
}

/* TODO store index of event ? */
/** catch variable events */
static
SCIP_RETCODE catchVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to catch bound change events */      
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdata->isremovedfixings = TRUE;
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->linbndchgeventdata, consdata->nlinvars) );
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      SCIPintervalSetEmpty(&consdata->linrange[i]);
      consdata->linbndchgeventdata[i].consdata = consdata;
      consdata->linbndchgeventdata[i].varidx = i;
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->linvars[i], SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr, &consdata->linbndchgeventdata[i], NULL) );
      consdata->isremovedfixings = consdata->isremovedfixings && SCIPvarIsActive(consdata->linvars[i]);
   }
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->quadbndchgeventdata, consdata->nquadvars) );
   SCIPintervalSetEmpty(&consdata->quadrange);
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      SCIPintervalSetEmpty(&consdata->quadrangevar[i]);
      consdata->quadbndchgeventdata[i].consdata = consdata;
      consdata->quadbndchgeventdata[i].varidx   = -i-1;
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->quadvars[i], SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr, &consdata->quadbndchgeventdata[i], NULL) );
      consdata->isremovedfixings = consdata->isremovedfixings && SCIPvarIsActive(consdata->quadvars[i]);
   }
   for( i = 0; i < consdata->nbilinterms; ++i )
      SCIPintervalSetEmpty(&consdata->bilinrange[i]);
   consdata->ispropagated = FALSE;

   return SCIP_OKAY;
}

/** drop variable events */
static
SCIP_RETCODE dropVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to drop bound change events */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( i = 0; i < consdata->nlinvars; ++i )
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->linvars[i],  SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr, &consdata->linbndchgeventdata[i], -1) );
   
   for( i = 0; i < consdata->nquadvars; ++i )
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->quadvars[i], SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr, &consdata->quadbndchgeventdata[i], -1) );

   return SCIP_OKAY;
}

/** sets or replaces function data of constraint
 * 
 *  Takes care of and release/capture of variables, but not of unlock/lock.
 */
static
SCIP_RETCODE consdataSetFunctionData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< the constraint data where to store function */
   SCIP_HASHMAP*         terms,              /**< linear and quadratic terms of function */
   int                   component           /**< which component of the terms, or -1 to use all */
   )
{  /*lint --e{777} */
   int                i, j;
   SCIP_HASHMAPLIST*  list;
   SCIP_HASHMAPLIST*  list2;
   SCIP_VAR*          var;
   SCIP_VAR*          bvar;
   PRESOLVEQUADTERM*  term;
   PRESOLVEBILINITEM* bitem;
   int                nlinvars  = 0;
   int                nquadvars = 0;
   int                nbilinterms   = 0;
   int                i_lin   = 0;
   int                i_quad  = 0;
   int                i_bilin = 0;
   int                i_adjbilin = 0;
   
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(terms != NULL);

   for( i = 0; i < consdata->nlinvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &consdata->linvars[i]) );
   }
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &consdata->quadvars[i]) );
   }

   /* collect statistics */
   for( i = 0; i < SCIPhashmapGetNLists(terms); ++i )
      for( list = SCIPhashmapGetList(terms, i); list != NULL; list = SCIPhashmapListGetNext(list) )
      {
         term = (PRESOLVEQUADTERM*) SCIPhashmapListGetImage(list);
         if( component >= 0 && term->component != component )
            continue;

         if( term->sqrcoefs == 0.0 && term->bilin == NULL )
         { /* linear variable */
            if( term->lincoefs != 0.0 )
               ++nlinvars;
            continue;
         }
         ++nquadvars;

         if( term->bilin != NULL )
            for( j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j )
               for( list2 = SCIPhashmapGetList(term->bilin, j); list2 != NULL; list2 = SCIPhashmapListGetNext(list2) )
                  ++nbilinterms;
      }

   nbilinterms /= 2; /* because we counted each bilinear term twice */

   /* realloc memory */
   consdata->nlinvars = nlinvars;
   if( nlinvars )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->linvars,   nlinvars) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->lincoefs, nlinvars) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->linrange, nlinvars) );
   }
   else
   {
      SCIPfreeMemoryArrayNull(scip, &consdata->linvars);
      SCIPfreeMemoryArrayNull(scip, &consdata->lincoefs);
      SCIPfreeMemoryArrayNull(scip, &consdata->linrange);
   }

   if( consdata->adjbilin != NULL )
      for( i = 0; i < consdata->nquadvars; ++i )
      {
         SCIPfreeMemoryArrayNull(scip, &consdata->adjbilin[i]);
      }

   consdata->nquadvars = nquadvars;
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->quadvars, nquadvars) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->quadlincoefs, nquadvars) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->quadsqrcoefs, nquadvars) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->nadjbilin, nquadvars) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->adjbilin, nquadvars) );

   consdata->nbilinterms = nbilinterms;
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->bilinvars1, nbilinterms) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->bilinvars2, nbilinterms) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->bilincoefs, nbilinterms) );

   /* set constraint data */
   for( i = 0; i < SCIPhashmapGetNLists(terms); ++i )
      for( list = SCIPhashmapGetList(terms, i); list != NULL; list = SCIPhashmapListGetNext(list) )
      {
         term = (PRESOLVEQUADTERM*) SCIPhashmapListGetImage(list);
         if( component >= 0 && term->component != component )
            continue;
         var  = (SCIP_VAR*) SCIPhashmapListGetOrigin(list);

         if( term->sqrcoefs == 0.0 && term->bilin == NULL )
         { /* linear variable */
            if( term->lincoefs != 0.0 )
            {
               assert(i_lin < nlinvars);
               assert(consdata->linvars   != NULL);
               assert(consdata->lincoefs != NULL);
               consdata->linvars  [i_lin] = var;
               consdata->lincoefs[i_lin] = term->lincoefs;
               ++i_lin;
            }
            continue;
         }

         /* quadratic variable */
         assert(i_quad < nquadvars);
         consdata->quadvars[i_quad] = var;
         consdata->quadlincoefs[i_quad] = term->lincoefs;
         consdata->quadsqrcoefs[i_quad] = term->sqrcoefs;
         consdata->nadjbilin[i_quad] = 0;
         if( term->bilin != NULL )
            for( j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j )
               for( list2 = SCIPhashmapGetList(term->bilin, j); list2 != NULL; list2 = SCIPhashmapListGetNext(list2) )
                  ++consdata->nadjbilin[i_quad];

         if( consdata->nadjbilin[i_quad] != 0 )
         { /* bilinear terms of quadratic variable */
            SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->adjbilin[i_quad], consdata->nadjbilin[i_quad]) );
            i_adjbilin = 0;
            for( j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j )
               for( list2 = SCIPhashmapGetList(term->bilin, j); list2 != NULL; list2 = SCIPhashmapListGetNext(list2) )
               {
                  bvar = (SCIP_VAR*) SCIPhashmapListGetOrigin(list2);
                  bitem = (PRESOLVEBILINITEM*) SCIPhashmapListGetImage(list2);

                  if( bitem->bilinidx >= 0 )
                  { /* bilinear term has been created already, just store its index in adjbilin */
                     consdata->adjbilin[i_quad][i_adjbilin] = bitem->bilinidx;
                     ++i_adjbilin;
                  }
                  else
                  { /* bilinear term needs to be created here */
                     PRESOLVEQUADTERM*  bterm;
                     PRESOLVEBILINITEM* bbitem;
                     bterm  = (PRESOLVEQUADTERM*) SCIPhashmapGetImage(terms, bvar);
                     assert(bterm != NULL);
                     assert(bterm->bilin != NULL);
                     bbitem = (PRESOLVEBILINITEM*) SCIPhashmapGetImage(bterm->bilin, var);
                     assert(bbitem != NULL);
                     assert(bbitem->bilinidx < 0);
                     assert(bbitem->coef == bitem->coef);

                     assert(i_bilin < nbilinterms);
                     consdata->bilinvars1 [i_bilin] = var;
                     consdata->bilinvars2 [i_bilin] = bvar;
                     consdata->bilincoefs[i_bilin] = bitem->coef;
                     consdata->adjbilin[i_quad][i_adjbilin] = i_bilin;

                     bitem->bilinidx = i_bilin;
                     bbitem->bilinidx = i_bilin;

                     ++i_adjbilin;
                     ++i_bilin;
                  }
               }
         }
         else
         {
            consdata->adjbilin[i_quad] = NULL;
         }

         ++i_quad;
      }

   assert(i_lin == nlinvars);
   assert(i_quad == nquadvars);
   assert(i_bilin == nbilinterms);

   SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->quadrangevar, nquadvars) );
   if( nbilinterms != 0 )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->bilinrange, nbilinterms) );
   }
   else
   {
      SCIPfreeMemoryArrayNull(scip, &consdata->bilinrange);
   }

   for( i = 0; i < consdata->nlinvars ; ++i )
   {
      SCIP_CALL( SCIPcaptureVar(scip, consdata->linvars[i]) );
   }
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      SCIP_CALL( SCIPcaptureVar(scip, consdata->quadvars[i]) );
   }

   return SCIP_OKAY;
}

/** frees constraint data structure */
static
void consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to constraint data to free */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   SCIPfreeMemoryArrayNull(scip, &(*consdata)->linvars);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->lincoefs);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->linrange);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->linbndchgeventdata);

   SCIPfreeMemoryArrayNull(scip, &(*consdata)->quadvars);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->quadlincoefs);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->quadsqrcoefs);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->quadrangevar);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->quadbndchgeventdata);

   if( (*consdata)->adjbilin != NULL )
   {
      int i;
      assert((*consdata)->nadjbilin != NULL);
      for( i = 0; i < (*consdata)->nquadvars; ++i )
      {
         SCIPfreeMemoryArrayNull(scip, &(*consdata)->adjbilin[i]);
      }
      SCIPfreeMemoryArray(scip, &(*consdata)->adjbilin);
      SCIPfreeMemoryArray(scip, &(*consdata)->nadjbilin);
   }

   SCIPfreeMemoryArrayNull(scip, &(*consdata)->bilinvars1);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->bilinvars2);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->bilincoefs);
   SCIPfreeMemoryArrayNull(scip, &(*consdata)->bilinrange);

   SCIPfreeMemory(scip, consdata);
   *consdata = NULL;
}

/** adds a variable with coefficients and bilinear term to the presolve data structure */
static
SCIP_RETCODE presolveQuadTermAdd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         terms,              /**< the terms where to add something */
   SCIP_VAR*             var,                /**< the variable to add */
   SCIP_Real             lincoefs,           /**< the linear coefficient of the variable */
   SCIP_Real             sqrcoefs,           /**< the coefficient of the square term of the variable */
   SCIP_VAR*             bilinvars,          /**< the other variable in a bilinear term, or NULL */  
   SCIP_Real             bilincoef           /**< the coefficient of the bilinear term */
   )
{
   PRESOLVEQUADTERM* term;
   
   assert(scip != NULL);
   assert(terms != NULL);
   assert(var != NULL);
   
   term = (PRESOLVEQUADTERM*) SCIPhashmapGetImage(terms, var);
   if( term == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &term) );
      term->lincoefs  = lincoefs;
      term->sqrcoefs  = sqrcoefs;
      term->bilin = NULL;
      term->component = -1;
      SCIP_CALL( SCIPhashmapInsert( terms, var, term ) );
   }
   else
   {
      term->lincoefs += lincoefs;
      term->sqrcoefs += sqrcoefs;
   }

   if( bilinvars != NULL )
   {
      PRESOLVEBILINITEM* bitem;
      if( term->bilin == NULL )
      {
         SCIP_CALL( SCIPhashmapCreate(&term->bilin, SCIPblkmem(scip), 10) ); /* TODO 10 a good value? */
         SCIP_CALL( SCIPallocBlockMemory(scip, &bitem) );
         SCIP_CALL( SCIPhashmapInsert(term->bilin, bilinvars, bitem) );
         bitem->coef = bilincoef;
         bitem->bilinidx = -1;
      }
      else
      {
         bitem = (PRESOLVEBILINITEM*) SCIPhashmapGetImage(term->bilin, bilinvars);
         if( bitem == NULL )
         {
            SCIP_CALL( SCIPallocBlockMemory(scip, &bitem) );
            SCIP_CALL( SCIPhashmapInsert(term->bilin, bilinvars, bitem) );
            bitem->coef = bilincoef;
            bitem->bilinidx = -1;
         }
         else
         {
            bitem->coef += bilincoef;
         }
      }
   }

   return SCIP_OKAY;
}

/** frees the presolve data */
static
void presolveQuadTermFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP**        terms               /**< quadratic terms to free */
   )
{
   SCIP_HASHMAPLIST*  list1;
   SCIP_HASHMAPLIST*  list2;
   PRESOLVEQUADTERM*  term;
   PRESOLVEBILINITEM* bitem;
   int i;
   int j;   

   assert(scip  != NULL);
   assert(terms != NULL);

   if( *terms == NULL )
      return;

   for( i = 0; i < SCIPhashmapGetNLists(*terms); ++i )
   {
      list1 = SCIPhashmapGetList(*terms, i);
      for( ; list1 != NULL; list1 = SCIPhashmapListGetNext(list1) )
      {
         term = (PRESOLVEQUADTERM*)SCIPhashmapListGetImage(list1);
         if( term->bilin != NULL )
         {
            for( j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j )
            {
               list2 = SCIPhashmapGetList(term->bilin, j);
               for( ; list2 != NULL; list2 = SCIPhashmapListGetNext(list2) )
               {
                  bitem = (PRESOLVEBILINITEM*)SCIPhashmapListGetImage(list2);
                  SCIPfreeBlockMemory(scip, &bitem);
               }
            }
            SCIPhashmapFree(&term->bilin);
         }
         SCIPfreeBlockMemory(scip, &term);
      }
   }
   SCIPhashmapFree(terms);
   assert(*terms == NULL);
}

/* Prints the quadratic function stored in the presolve data structure. */
#ifdef SCIP_DEBUG
static
void presolveQuadTermPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         terms,              /**< quadratic terms */
   FILE*                 file                /**< file to print to, or NULL for stdout */
   )
{
   SCIP_HASHMAPLIST*  list1;
   SCIP_HASHMAPLIST*  list2;
   SCIP_VAR*          var;
   SCIP_VAR*          bvar;
   PRESOLVEQUADTERM*  term;
   PRESOLVEBILINITEM* bitem;
   int                i, j;
   
   assert(scip != NULL);
   assert(terms != NULL);

   for( i = 0; i < SCIPhashmapGetNLists(terms); ++i )
   {
      list1 = SCIPhashmapGetList(terms, i);
      for( ; list1; list1 = SCIPhashmapListGetNext(list1) )
      {
         var = (SCIP_VAR*)SCIPhashmapListGetOrigin(list1);
         term = (PRESOLVEQUADTERM*)SCIPhashmapListGetImage (list1);
         if( term->lincoefs != 0.0 )
            SCIPinfoMessage(scip, file, "%+g*%s ",      term->lincoefs, SCIPvarGetName(var));
         if( term->sqrcoefs != 0.0 )
            SCIPinfoMessage(scip, file, "%+g*sqr(%s) ", term->sqrcoefs, SCIPvarGetName(var));
         if( term->bilin != NULL )
         {
            SCIPinfoMessage(scip, file, "+%s*(", SCIPvarGetName(var));
            for( j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j )
            {
               list2 = SCIPhashmapGetList(term->bilin, j);
               for( ; list2 != NULL; list2 = SCIPhashmapListGetNext(list2) )
               {
                  bvar = (SCIP_VAR*)SCIPhashmapListGetOrigin(list2);
                  bitem = (PRESOLVEBILINITEM*)SCIPhashmapListGetImage (list2);
                  if( var < bvar )
                     SCIPinfoMessage(scip, file, "%+g*%s ", bitem->coef, SCIPvarGetName(bvar));
               }
            }
            SCIPinfoMessage(scip, file, ")");
         }
      }
   }
}
#endif

/** adds a linear term into the presolve data structure. */
static
SCIP_RETCODE presolveAddLinearTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         terms,              /**< presolve quadratic terms */
   SCIP_Real*            constant,           /**< presolve constant term */
   SCIP_Real             coef,               /**< coefficient of variable */
   SCIP_VAR*             var                 /**< variable to add */
   )
{
   assert(scip != NULL);
   assert(terms != NULL);
   assert(constant != NULL);
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);

   if( coef == 0.0 )
      return SCIP_OKAY;

   switch( SCIPvarGetStatus(var) )
   {
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
      {
         if( SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) ) /* var is fixed */
            *constant += coef * SCIPvarGetLbGlobal(var);
         else
         {
            SCIP_CALL( presolveQuadTermAdd(scip, terms, var, coef, 0.0, NULL, 0.0) );
         }
         break;
      }
      case SCIP_VARSTATUS_FIXED:
      {
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)));
         *constant += coef * SCIPvarGetLbGlobal(var);
         break;
      }
      case SCIP_VARSTATUS_AGGREGATED:
      { /* var is replaced by scalar * aggrvar + constant */
         *constant += coef * SCIPvarGetAggrConstant(var);
         SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, coef * SCIPvarGetAggrScalar(var), SCIPvarGetAggrVar(var)) );
         break;
      }
      case SCIP_VARSTATUS_MULTAGGR:
      { /* var is replaced by sum_i scalar_i * aggrvar_i + constant */
         int navars;
         SCIP_VAR** avars;
         SCIP_Real* ascalars;
         int j;

         navars = SCIPvarGetMultaggrNVars(var);
         avars = SCIPvarGetMultaggrVars(var);
         ascalars = SCIPvarGetMultaggrScalars(var);

         SCIPdebugMessage("replace multiaggregated variable %s\n", SCIPvarGetName(var));
         *constant += coef * SCIPvarGetMultaggrConstant(var);
         for( j = 0; j < navars; ++j )
         {
            SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, coef * ascalars[j], avars[j]) );
         }
         break;
      }
      case SCIP_VARSTATUS_NEGATED:
      { /* var is replaced by constant - negvar */
         *constant += coef * SCIPvarGetNegationConstant(var);
         SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, -coef, SCIPvarGetNegationVar(var)) );
         break;
      }
      case SCIP_VARSTATUS_ORIGINAL: /* for lint */
      default:
         SCIPerrorMessage("unexpected variable status: %d\n", SCIPvarGetStatus(var));
         return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** adds a bilinear term into the presolve data structure */  
static
SCIP_RETCODE presolveAddBilinearTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         terms,              /**< presolve quadratic terms */
   SCIP_Real*            constant,           /**< presolve constant term */
   SCIP_Real             coef,               /**< coefficient of bilinear variable */
   SCIP_VAR*             var1,               /**< first variable in bilinear term to add */
   SCIP_VAR*             var2                /**< second variable in bilinear term to add */
   )
{
   assert(scip != NULL);
   assert(terms != NULL);
   assert(constant != NULL);
   assert(var1 != NULL);
   assert(var2 != NULL);
   assert(SCIPvarGetStatus(var1) != SCIP_VARSTATUS_ORIGINAL);
   assert(SCIPvarGetStatus(var2) != SCIP_VARSTATUS_ORIGINAL);

   if( coef == 0.0 )
      return SCIP_OKAY;

   switch( SCIPvarGetStatus(var2) )
   {
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
      {
         if( SCIPisEQ(scip, SCIPvarGetLbGlobal(var2), SCIPvarGetUbGlobal(var2)) )
         { /* var2 is fixed */
            SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, coef*SCIPvarGetLbGlobal(var2), var1) );
         }
         else if( var1 == var2 )
         { /* var2 is not fixed but equal to var1 */
            SCIP_CALL( presolveQuadTermAdd(scip, terms, var1, 0.0, coef, NULL, 0.0) );
         }
         else
         { /* var2 is not fixed */
            assert(var1 != var2);
            if( SCIPvarIsActive(var1) )
            {
               SCIP_CALL( presolveQuadTermAdd(scip, terms, var1, 0.0, 0.0, var2, coef) );
               SCIP_CALL( presolveQuadTermAdd(scip, terms, var2, 0.0, 0.0, var1, coef) );
            }
            else
            {
               SCIP_CALL( presolveAddBilinearTerm(scip, terms, constant, coef, var2, var1) );
            }
         }
         break;
      }
      case SCIP_VARSTATUS_FIXED:
      {
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var2), SCIPvarGetUbGlobal(var2)));
         SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, coef*SCIPvarGetLbGlobal(var2), var1) );
         break;
      }
      case SCIP_VARSTATUS_AGGREGATED:
      { /* var2 is replaced by scalar * aggrvar + constant */
         SCIP_CALL( presolveAddLinearTerm  (scip, terms, constant, coef*SCIPvarGetAggrConstant(var2), var1) );
         SCIP_CALL( presolveAddBilinearTerm(scip, terms, constant, coef*SCIPvarGetAggrScalar(var2),   var1, SCIPvarGetAggrVar(var2)) );
         break;
      }
      case SCIP_VARSTATUS_MULTAGGR:
      { /* var2 is replaced by sum_i scalar_i * aggrvar_i + constant */
         int navars;
         SCIP_VAR** avars;
         SCIP_Real* ascalars;
         int j;

         navars   = SCIPvarGetMultaggrNVars(var2);
         avars    = SCIPvarGetMultaggrVars(var2);
         ascalars = SCIPvarGetMultaggrScalars(var2);

         SCIPdebugMessage("replace multiaggregated variable %s in term %g*%s*%s by %g", 
            SCIPvarGetName(var2), coef, SCIPvarGetName(var1), SCIPvarGetName(var2), SCIPvarGetMultaggrConstant(var2));

         for( j = 0; j < navars; ++j )
         {
            SCIPdebugPrintf(" + %g%s", ascalars[j], SCIPvarGetName(avars[j]));
            SCIP_CALL( presolveAddBilinearTerm(scip, terms, constant, coef * ascalars[j], var1, avars[j]) );
         }
         SCIPdebugPrintf("\n");
         SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, coef * SCIPvarGetMultaggrConstant(var2), var1) );
         break;
      }
      case SCIP_VARSTATUS_NEGATED:
      { /* var2 is replaced by constant - negvar */
         SCIP_CALL( presolveAddLinearTerm  (scip, terms, constant, coef * SCIPvarGetNegationConstant(var2), var1) );
         SCIP_CALL( presolveAddBilinearTerm(scip, terms, constant, -coef, var1, SCIPvarGetNegationVar(var2)) );
         break;
      }
      case SCIP_VARSTATUS_ORIGINAL: /* for lint */
      default:
         SCIPerrorMessage("unexpected variable status: %d\n", SCIPvarGetStatus(var2));
         return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** adds a square term into the presolve data structure. */
static
SCIP_RETCODE presolveAddSquareTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         terms,              /**< presolve quadratic terms */
   SCIP_Real*            constant,           /**< presolve constant term */
   SCIP_Real             coef,               /**< coefficient of variable */
   SCIP_VAR*             var                 /**< variable to add */
   )
{
   assert(scip != NULL);
   assert(terms != NULL);
   assert(constant != NULL);
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);

   if( coef == 0.0 )
      return SCIP_OKAY;

   switch( SCIPvarGetStatus(var) )
   {
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_COLUMN:
      {
         if( SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) )
         { /* variable is fixed now */
            SCIP_Real val;

            val = SCIPvarGetLbGlobal(var);
            *constant += coef * val * val;
         }
         else
         { /* variable is not fixed */
            SCIP_CALL( presolveQuadTermAdd(scip, terms, var, 0.0, coef, NULL, 0.0) );
         }
         break;
      }
      case SCIP_VARSTATUS_FIXED:
      {
         SCIP_Real val;
         val = SCIPvarGetLbGlobal(var);
         assert(SCIPisEQ(scip, val, SCIPvarGetUbGlobal(var)));
         *constant += coef * val * val;
         break;
      }
      case SCIP_VARSTATUS_AGGREGATED:
      { /* var is replaced by scalar * aggrvar + constant */
         SCIP_Real aggconst;
         SCIP_Real aggscalar;

         aggconst = SCIPvarGetAggrConstant(var);
         aggscalar = SCIPvarGetAggrScalar(var);

         SCIP_CALL( presolveAddSquareTerm(scip, terms, constant,   coef*aggscalar*aggscalar, SCIPvarGetAggrVar(var)) );
         SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, 2*coef*aggscalar*aggconst , SCIPvarGetAggrVar(var)) );
         *constant += coef * aggconst * aggconst;
         break;
      }
      case SCIP_VARSTATUS_MULTAGGR:
      { /* var is replaced by sum_i scalar_i * aggrvar_i + constant */
         SCIP_VAR** avars;
         SCIP_Real* ascalars;
         SCIP_Real aconstant;
         int navars;
         int j;
         int k;

         navars = SCIPvarGetMultaggrNVars(var);
         avars = SCIPvarGetMultaggrVars(var);
         ascalars = SCIPvarGetMultaggrScalars(var);
         aconstant = SCIPvarGetMultaggrConstant(var);

         SCIPdebugMessage("replace multiaggregated variable %s by %g", SCIPvarGetName(var), aconstant);
         *constant += coef * aconstant * aconstant;
         for( j = 0; j < navars; ++j )
         {
            SCIPdebugPrintf("+ %g%s", ascalars[j], SCIPvarGetName(avars[j]));
            SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, 2*coef*ascalars[j]*aconstant, avars[j]) );
            for( k = 0; k < j; ++k )
            {
               SCIP_CALL( presolveAddBilinearTerm(scip, terms, constant, 2*coef*ascalars[j]*ascalars[k], avars[k], avars[j]) );
            }
            SCIP_CALL( presolveAddSquareTerm(scip, terms, constant, coef*ascalars[j]*ascalars[j], avars[j]) );
         }
         SCIPdebugPrintf("\n");
         break;
      }
      case SCIP_VARSTATUS_NEGATED:
      { /* var is replaced by constant - negvar */
         SCIP_Real negconst;

         negconst = SCIPvarGetNegationConstant(var);
         *constant += coef * negconst * negconst;
         SCIP_CALL( presolveAddLinearTerm(scip, terms, constant, -2*coef*negconst, SCIPvarGetNegationVar(var)) );
         SCIP_CALL( presolveAddSquareTerm(scip, terms, constant, coef, SCIPvarGetNegationVar(var)) );
         break;
      }
      case SCIP_VARSTATUS_ORIGINAL: /* for lint */
      default:
         SCIPerrorMessage("unexpected variable status: %d\n", SCIPvarGetStatus(var));
         return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** gets constraint function data as set of PRESOLVEQUADTERM's for reformulation and preprocessing
 * 
 *  replaces fixed and aggregated variables;
 *  replaces squares of binary variables
 */
static
SCIP_RETCODE presolveCreateQuadTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP**        terms,              /**< storage for PRESOLVEQUADTERM's */
   SCIP_Real*            constant,           /**< storage for constant */
   SCIP_Bool*            havechange,         /**< was a change in the function has been identified (variables fixed, aggregated,...) ? */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             replacesqrbinary    /**< should we replace squares of binary variables ? */
   )
{
   SCIP_CONSDATA*     consdata;
   SCIP_VAR*          var;
   SCIP_Real          coef;
   SCIP_HASHMAPLIST*  listitem;
   SCIP_HASHMAPLIST*  blistitem;
   PRESOLVEQUADTERM*  term;
   PRESOLVEBILINITEM* bilinitem;
   SCIP_Bool          havebilin;
   int                i;
   int                j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(terms != NULL);
   assert(constant != NULL);
   assert(havechange != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( SCIPhashmapCreate(terms, SCIPblkmem(scip), consdata->nquadvars + consdata->nlinvars) );
   *constant = 0.;
   *havechange = FALSE;

   for( i = 0; i < consdata->nlinvars; ++i )
   {
      var = consdata->linvars[i];
      *havechange = *havechange || ((SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN) && (SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE));
      *havechange = *havechange || SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
      *havechange = *havechange || SCIPisZero(scip, consdata->lincoefs[i]);
      SCIP_CALL( presolveAddLinearTerm(scip, *terms, constant, consdata->lincoefs[i], var) );
   }

   for( j = 0; j < consdata->nquadvars; ++j )
   {
      if( consdata->quadlincoefs[j] == 0.0 && consdata->quadsqrcoefs[j] == 0.0 )
      {
         if( consdata->nadjbilin[j] == 0 )
            *havechange = TRUE;
         continue;
      }

      var = consdata->quadvars[j];
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL); /* should not happen after transformation */
      switch( SCIPvarGetStatus(var) )
      {
         case SCIP_VARSTATUS_LOOSE:
         case SCIP_VARSTATUS_COLUMN:
         {
            if( SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) )
            { /* variable is fixed now */
               SCIP_Real val;

               val = SCIPvarGetLbGlobal(var);
               *constant += consdata->quadlincoefs[j] * val;
               *constant += consdata->quadsqrcoefs[j] * val*val;
               *havechange = TRUE;
            }
            else
            { /* variable is not fixed */
               SCIP_CALL( presolveQuadTermAdd(scip, *terms, var, consdata->quadlincoefs[j], consdata->quadsqrcoefs[j], NULL, 0.0) );
            }
            break;
         }
         case SCIP_VARSTATUS_FIXED:
         {
            SCIP_Real val;
            
            val = SCIPvarGetLbGlobal(var);
            assert(SCIPisEQ(scip, val, SCIPvarGetUbGlobal(var)));
            *constant += consdata->quadlincoefs[j] * val;
            *constant += consdata->quadsqrcoefs[j] * val*val;
            *havechange = TRUE;
            break;
         }
         case SCIP_VARSTATUS_AGGREGATED:
         case SCIP_VARSTATUS_MULTAGGR:
         case SCIP_VARSTATUS_NEGATED:
         {
            SCIP_CALL( presolveAddLinearTerm(scip, *terms, constant, consdata->quadlincoefs[j], var) );
            SCIP_CALL( presolveAddSquareTerm(scip, *terms, constant, consdata->quadsqrcoefs[j], var) );
            *havechange = TRUE;
            break;
         }
         case SCIP_VARSTATUS_ORIGINAL: /* for lint */
         default:
            SCIPerrorMessage("unexpected variable status: %d\n", SCIPvarGetStatus(var));
            return SCIP_ERROR;
      }
   }

   for( j = 0; j < consdata->nbilinterms; ++j )
   { /* go through bilinear terms */
      var   = consdata->bilinvars1[j];
      coef = consdata->bilincoefs[j];

      if( coef == 0.0 )
      {
         *havechange = TRUE;
         continue;
      }

      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL);
      switch( SCIPvarGetStatus(var) )
      {
         case SCIP_VARSTATUS_LOOSE:
         case SCIP_VARSTATUS_COLUMN:
         {
            if( SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) )
            { /* variable is now fixed */
               SCIP_CALL( presolveAddLinearTerm(scip, *terms, constant, coef*SCIPvarGetLbGlobal(var), consdata->bilinvars2[j]) );
               *havechange = TRUE;
            }
            else
            { /* variable is not fixed */
               *havechange = *havechange || ((SCIPvarGetStatus(consdata->bilinvars2[j]) != SCIP_VARSTATUS_COLUMN) && (SCIPvarGetStatus(consdata->bilinvars2[j]) != SCIP_VARSTATUS_LOOSE));
               *havechange = *havechange || SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->bilinvars2[j]), SCIPvarGetUbGlobal(consdata->bilinvars2[j]));
               SCIP_CALL( presolveAddBilinearTerm(scip, *terms, constant, coef, var, consdata->bilinvars2[j]) );
            }
            break;
         }
         case SCIP_VARSTATUS_FIXED:
         case SCIP_VARSTATUS_AGGREGATED:
         case SCIP_VARSTATUS_MULTAGGR:
         case SCIP_VARSTATUS_NEGATED:
         {
            SCIP_CALL( presolveAddBilinearTerm(scip, *terms, constant, consdata->bilincoefs[j], consdata->bilinvars2[j], var) );
            *havechange = TRUE;
            break;
         }
         case SCIP_VARSTATUS_ORIGINAL: /* for lint */
         default:
            SCIPerrorMessage("unexpected variable status: %d\n", SCIPvarGetStatus(var));
            return SCIP_ERROR;
      }
   }

   /* clean 0's */
   for( i = 0; i < SCIPhashmapGetNLists(*terms); ++i )
   {
      listitem = SCIPhashmapGetList(*terms, i);
      while( listitem != NULL )
      {
         var  = (SCIP_VAR*) SCIPhashmapListGetOrigin(listitem);
         term = (PRESOLVEQUADTERM*) SCIPhashmapListGetImage(listitem);
         if( replacesqrbinary && SCIPvarGetType(var) == SCIP_VARTYPE_BINARY && term->sqrcoefs != 0.0 )
         { /* replace square of binary variable by variable itself (x^2 = x) */
            term->lincoefs += term->sqrcoefs;
            term->sqrcoefs  = 0.0;
            *havechange = TRUE;
         }

         if( SCIPisZero(scip, term->lincoefs) )
            term->lincoefs = 0.0;
         if( SCIPisZero(scip, term->sqrcoefs) )
            term->sqrcoefs = 0.0;

         havebilin = FALSE;
         if( term->bilin )
         {
            for( j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j )
            {
               blistitem = SCIPhashmapGetList(term->bilin, j);
               while( blistitem != NULL )
               {
                  bilinitem = (PRESOLVEBILINITEM*) SCIPhashmapListGetImage(blistitem);
                  if( SCIPisZero(scip, bilinitem->coef) )
                  {
                     SCIPfreeBlockMemory(scip, &bilinitem);
                     SCIP_CALL( SCIPhashmapRemove(term->bilin, SCIPhashmapListGetOrigin(blistitem)) );
                     blistitem = SCIPhashmapGetList(term->bilin, j); /* restart from the beginning of this list again */
                  }
                  else
                  {
                     blistitem  = SCIPhashmapListGetNext(blistitem);
                     havebilin = TRUE;
                  }
               }
            }
            if( !havebilin )
               SCIPhashmapFree(&term->bilin);
         }

         if( term->lincoefs == 0.0 && term->sqrcoefs == 0.0 && !havebilin)
         {
            SCIPfreeBlockMemory(scip, &term);
            SCIP_CALL( SCIPhashmapRemove(*terms, SCIPhashmapListGetOrigin(listitem)) );
            listitem = SCIPhashmapGetList(*terms, i); /* restart from the beginning of this list again */
         }
         else
            listitem = SCIPhashmapListGetNext(listitem);
      }
   }

   return SCIP_OKAY;
}

/** assigns component (block) number to all variables associated to a given one */
static
void presolveAssignComponent(
   SCIP_HASHMAP*         terms,              /**< PRESOLVEQUADTERM's */
   SCIP_VAR*             var,                /**< variable to assign to a component */
   int                   component           /**< component number to assign */
   )
{
   PRESOLVEQUADTERM* term;
   SCIP_HASHMAPLIST* list;
   int               j;
   
   assert(terms != NULL);
   assert(var != NULL);
   
   term = (PRESOLVEQUADTERM*)SCIPhashmapGetImage(terms, var);
   assert(term != NULL);
   assert(term->component == component || term->component == -1);

   if( term->component >= 0 )
      return;
   
   term->component = component;

   if( term->bilin == NULL )
      return;
   
   for( j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j )
      for( list = SCIPhashmapGetList(term->bilin, j); list != NULL; list = SCIPhashmapListGetNext(list) )
         presolveAssignComponent(terms, (SCIP_VAR*)SCIPhashmapListGetOrigin(list), component);
}

/** find components (blocks) of connected quadratic variables in constraint function data
 */
static
SCIP_RETCODE presolveFindComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         terms,              /**< constraint function in form of PRESOLVEQUADTERM's */
   int*                  ncomponents         /**< buffer where to store number of components found */
   )
{
   SCIP_HASHMAPLIST* list;
   SCIP_HASHMAPLIST* list2;
   PRESOLVEQUADTERM* term;
   int i;
   int j;
   
   assert(scip != NULL);
   assert(terms != NULL);
   assert(ncomponents != NULL);

   *ncomponents = 0;
   for( i = 0; i < SCIPhashmapGetNLists(terms); ++i )
      for( list = SCIPhashmapGetList(terms, i); list != NULL; list = SCIPhashmapListGetNext(list) )
      {
         term = (PRESOLVEQUADTERM*) SCIPhashmapListGetImage(list);
         if( term->sqrcoefs == 0.0 && term->bilin == NULL ) /* linear variable */
            continue;
         if( term->component >= 0 ) /* component already assigned */
            continue;
         
         term->component = *ncomponents;
         ++*ncomponents;
         
         if( term->bilin == NULL )
            continue;
         for( j = 0; j < SCIPhashmapGetNLists(term->bilin); ++j )
            for( list2 = SCIPhashmapGetList(term->bilin, j); list2 != NULL; list2 = SCIPhashmapListGetNext(list2) )
               presolveAssignComponent(terms, (SCIP_VAR*)SCIPhashmapListGetOrigin(list2), term->component);
      }
   
   return SCIP_OKAY;
}

#if 1
/** disaggregates a block separable constraint into several quadratic constraints
 */
static
SCIP_RETCODE presolveDisaggregate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler itself */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_HASHMAP*         terms,              /**< constraint function in form of PRESOLVEQUADTERM's */
   SCIP_Real             constant,           /**< constant part of constraint function */
   int                   ncomponents         /**< number of components (blocks) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_HASHMAPLIST*  list;
   PRESOLVEQUADTERM*  term;
   SCIP_VAR**         newlinvars;
   SCIP_Real*         newlincoefs;
   SCIP_CONS*         blockcons;
   SCIP_CONSDATA*     blockconsdata;
   SCIP_VAR*          auxvar;
   SCIP_VAR*          var;
   int                nnewlinvars;
   int                inewlinvars;
   char               name[255];
   int                i;
   int                k;
   SCIP_CONS*         lincons;
   SCIP_Real          lhs, rhs;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(terms != NULL);
   assert(ncomponents > 1);
   assert(SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("disaggregate constraint %s\n", SCIPconsGetName(cons));

   nnewlinvars = ncomponents;
   for( i = 0; i < SCIPhashmapGetNLists(terms); ++i )
      for( list = SCIPhashmapGetList(terms, i); list != NULL; list = SCIPhashmapListGetNext(list) )
      {
         term = (PRESOLVEQUADTERM*) SCIPhashmapListGetImage(list);
         if( term->component == -1 )
         {
            assert(SCIPisZero(scip, term->sqrcoefs));
            assert(term->bilin == NULL);
            if( SCIPisZero(scip, term->lincoefs) )
               continue;
            ++nnewlinvars;
         }
      }

   SCIP_CALL( SCIPallocBufferArray(scip, &newlinvars,  nnewlinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newlincoefs, nnewlinvars) );

   inewlinvars = 0;
   for( k = 0; k < ncomponents; ++k )
   {
      SCIP_CALL( SCIPallocMemory( scip, &blockconsdata) );

      /* we need to enforce only one bound here
       * we cannot make auxvar implicit integer then
       * hope it does not have bad effect on bound tightening
       */
      blockconsdata->lhs = SCIPisInfinity(scip, -consdata->lhs) ? -SCIPinfinity(scip) : 0.;
      blockconsdata->rhs = SCIPisInfinity(scip,  consdata->rhs) ?  SCIPinfinity(scip) : 0.;

      blockconsdata->nlinvars = 0;
      blockconsdata->linvars = NULL;
      blockconsdata->lincoefs = NULL;

      blockconsdata->nquadvars = 0;
      blockconsdata->quadvars = NULL;
      blockconsdata->quadlincoefs = NULL;
      blockconsdata->quadsqrcoefs = NULL;
      blockconsdata->nadjbilin = NULL;
      blockconsdata->adjbilin = NULL;

      blockconsdata->nbilinterms = 0;
      blockconsdata->bilincoefs = NULL;
      blockconsdata->bilinvars1 = NULL;
      blockconsdata->bilinvars2 = NULL;

      blockconsdata->linrange = NULL;
      blockconsdata->quadrangevar = NULL;
      blockconsdata->linbndchgeventdata = NULL;
      blockconsdata->quadbndchgeventdata = NULL;
      blockconsdata->bilinrange = NULL;

      SCIP_CALL( consdataSetFunctionData(scip, blockconsdata, terms, k) );
      assert(blockconsdata->nlinvars == 0);
      assert(blockconsdata->linvars == NULL);
      assert(blockconsdata->lincoefs == NULL);
      assert(blockconsdata->linrange == NULL);
      /* we scale the new constraints by ncomponents, so that the sum of the feasibility violations stays below the feasibility tolerance in a solution */
#if 0
      for( i = 0; i < blockconsdata->nquadvars; ++i )
      {
         blockconsdata->quadlincoefs[i] *= ncomponents;
         blockconsdata->quadsqrcoefs[i] *= ncomponents;
      }
      for( i = 0; i < blockconsdata->nbilinterms; ++i )
         blockconsdata->bilincoefs[i] *= ncomponents;
#endif

      (void)SCIPsnprintf(name, 255, "%s#%u", SCIPconsGetName(cons), k);
      SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS, SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvar) );

      SCIP_CALL( SCIPallocMemoryArray(scip, &blockconsdata->linvars,  1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &blockconsdata->lincoefs, 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &blockconsdata->linrange, 1) );

      blockconsdata->nlinvars = 1;
      blockconsdata->linvars[0] = auxvar;
      blockconsdata->lincoefs[0] = -1;

      blockconsdata->isconvex = FALSE;
      blockconsdata->isconcave = FALSE;
      blockconsdata->isremovedfixings = TRUE;
      blockconsdata->ispropagated = FALSE;
      blockconsdata->ispresolved = FALSE; /* so that in the next presolve round maybe auxvar is made implicit integer */

      SCIP_CALL( SCIPcreateCons(scip, &blockcons, name, conshdlr, blockconsdata, 
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
         SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
         SCIPconsIsStickingAtNode(cons)) );

      SCIP_CALL( SCIPaddCons(scip, blockcons) );
      SCIPdebugMessage("created new constraint %s: ", SCIPconsGetName(blockcons));
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintCons(scip, blockcons, NULL) );
#endif
      if( SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVING )
      {
         SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, blockcons) );
      }

      SCIP_CALL( SCIPreleaseCons(scip, &blockcons) );

      assert(inewlinvars < nnewlinvars);
      newlinvars [inewlinvars] = auxvar;
#if 0
      newlincoefs[inewlinvars] = 1.0/ncomponents;
#else
      newlincoefs[inewlinvars] = 1.0;
#endif
      ++inewlinvars;
   }

   for( i = 0; i < SCIPhashmapGetNLists(terms); ++i )
      for( list = SCIPhashmapGetList(terms, i); list != NULL; list = SCIPhashmapListGetNext(list) )
      {
         term = (PRESOLVEQUADTERM*) SCIPhashmapListGetImage(list);
         if( term->component == -1 && !SCIPisZero(scip, term->lincoefs) )
         {
            var = (SCIP_VAR*) SCIPhashmapListGetOrigin(list);
            assert(inewlinvars < nnewlinvars);
            newlinvars [inewlinvars] = var;
            newlincoefs[inewlinvars] = term->lincoefs;
            ++inewlinvars;
         }
      }
   assert(inewlinvars == nnewlinvars);

   /* create linear constraint from linear part and auxiliary variables */
   if( !SCIPisInfinity(scip, -consdata->lhs) )
      lhs = consdata->lhs - constant;
   else
      lhs = -SCIPinfinity(scip);
   if( !SCIPisInfinity(scip, consdata->rhs) )
      rhs = consdata->rhs - constant;
   else
      rhs = SCIPinfinity(scip);
   
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons), nnewlinvars, newlinvars, newlincoefs, lhs, rhs,
      SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
      SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
      SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
      SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   
   SCIPdebugMessage("added coupling constraint ");
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, lincons, NULL) ) );
   
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

   SCIPfreeBufferArray(scip, &newlinvars);
   SCIPfreeBufferArray(scip, &newlincoefs);

   return SCIP_OKAY;
}
#else
/** disaggregates a block separable constraint into several quadratic constraints
 */
static
SCIP_RETCODE presolveDisaggregate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler itself */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_HASHMAP*         terms,              /**< constraint function in form of PRESOLVEQUADTERM's */
   SCIP_Real             constant,           /**< constant part of constraint function */
   int                   ncomponents         /**< number of components (blocks) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_HASHMAPLIST*  list;
   PRESOLVEQUADTERM*  term;
   SCIP_VAR**         newlinvars;
   SCIP_Real*         newlincoefs;
   SCIP_CONS*         blockcons;
   SCIP_CONSDATA*     blockconsdata;
   SCIP_VAR*          auxvar;
   SCIP_VAR*          var;
   int                nnewlinvars;
   int                inewlinvars;
   char               name[255];
   int                i;
   int                k;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(terms != NULL);
   assert(ncomponents > 1);
   assert(SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("disaggregate constraint %s\n", SCIPconsGetName(cons));

   nnewlinvars = ncomponents - 1;
   for( i = 0; i < SCIPhashmapGetNLists(terms); ++i )
      for( list = SCIPhashmapGetList(terms, i); list != NULL; list = SCIPhashmapListGetNext(list) )
      {
         term = (PRESOLVEQUADTERM*) SCIPhashmapListGetImage(list);
         if( term->component == -1 )
         {
            assert(SCIPisZero(scip, term->sqrcoefs));
            assert(term->bilin == NULL);
            if( SCIPisZero(scip, term->lincoefs) )
               continue;
            ++nnewlinvars;
         }
      }

   SCIP_CALL( SCIPallocMemoryArray(scip, &newlinvars,   nnewlinvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &newlincoefs, nnewlinvars) );

   inewlinvars = 0;
   for( k = 1; k < ncomponents; ++k )
   {
      SCIP_CALL( SCIPallocMemory( scip, &blockconsdata) );

      /* we need to enforce only one bound here
       * we cannot make auxvar implicit integer then
       * hope it does not have bad effect on bound tightening
       */
      blockconsdata->lhs = SCIPisInfinity(scip, -consdata->lhs) ? -SCIPinfinity(scip) : 0.;
      blockconsdata->rhs = SCIPisInfinity(scip,  consdata->rhs) ?  SCIPinfinity(scip) : 0.;

      blockconsdata->nlinvars = 0;
      blockconsdata->linvars = NULL;
      blockconsdata->lincoefs = NULL;

      blockconsdata->nquadvars = 0;
      blockconsdata->quadvars = NULL;
      blockconsdata->quadlincoefs = NULL;
      blockconsdata->quadsqrcoefs = NULL;
      blockconsdata->nadjbilin = NULL;
      blockconsdata->adjbilin = NULL;

      blockconsdata->nbilinterms = 0;
      blockconsdata->bilincoefs = NULL;
      blockconsdata->bilinvars1 = NULL;
      blockconsdata->bilinvars2 = NULL;

      blockconsdata->linrange = NULL;
      blockconsdata->quadrangevar = NULL;
      blockconsdata->linbndchgeventdata = NULL;
      blockconsdata->quadbndchgeventdata = NULL;
      blockconsdata->bilinrange = NULL;

      SCIP_CALL( consdataSetFunctionData(scip, blockconsdata, terms, k) );
      assert(blockconsdata->nlinvars == 0);
      assert(blockconsdata->linvars == NULL);
      assert(blockconsdata->lincoefs == NULL);
      assert(blockconsdata->linrange == NULL);
      /* we scale the new constraints by ncomponents, so that the sum of the feasibility violations stays below the feasibility tolerance in a solution */
#if 0
      for( i = 0; i < blockconsdata->nquadvars; ++i )
      {
         blockconsdata->quadlincoefs[i] *= ncomponents;
         blockconsdata->quadsqrcoefs[i] *= ncomponents;
      }
      for( i = 0; i < blockconsdata->nbilinterms; ++i )
         blockconsdata->bilincoefs[i] *= ncomponents;
#endif

      (void)SCIPsnprintf(name, 255, "%s#%u", SCIPconsGetName(cons), k);
      SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS, SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvar) );

      SCIP_CALL( SCIPallocMemoryArray(scip, &blockconsdata->linvars,   1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &blockconsdata->lincoefs, 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &blockconsdata->linrange, 1) );

      blockconsdata->nlinvars = 1;
      blockconsdata->linvars[0] = auxvar;
      blockconsdata->lincoefs[0] = -1;

      blockconsdata->isconvex = FALSE;
      blockconsdata->isconcave = FALSE;
      blockconsdata->isremovedfixings = TRUE;
      blockconsdata->ispropagated = FALSE;
      blockconsdata->ispresolved = FALSE; /* so that in the next presolve round maybe auxvar is made implicit integer */

      SCIP_CALL( SCIPcreateCons(scip, &blockcons, name, conshdlr, blockconsdata, 
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
         SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
         SCIPconsIsStickingAtNode(cons)) );

      SCIP_CALL( SCIPaddCons(scip, blockcons) );
      SCIPdebugMessage("created new constraint %s: ", SCIPconsGetName(blockcons));
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintCons(scip, blockcons, NULL) );
#endif
      if( SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVING )
      {
         SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, blockcons) );
      }

      SCIP_CALL( SCIPreleaseCons(scip, &blockcons) );

      assert(inewlinvars < nnewlinvars);
      newlinvars  [inewlinvars] = auxvar;
#if 0
      newlincoefs[inewlinvars] = 1.0/ncomponents;
#else
      newlincoefs[inewlinvars] = 1.0;
#endif
      ++inewlinvars;
      SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   }

   SCIP_CALL( consdataSetFunctionData(scip, consdata, terms, 0) );
   assert(consdata->nlinvars == 0);
   assert(consdata->linvars == NULL);
   assert(consdata->lincoefs == NULL);
   assert(consdata->linrange == NULL);

   for( i = 0; i < SCIPhashmapGetNLists(terms); ++i )
      for( list = SCIPhashmapGetList(terms, i); list != NULL; list = SCIPhashmapListGetNext(list) )
      {
         term = (PRESOLVEQUADTERM*) SCIPhashmapListGetImage(list);
         if( term->component == -1 && !SCIPisZero(scip, term->lincoefs) )
         {
            var = (SCIP_VAR*) SCIPhashmapListGetOrigin(list);
            assert(inewlinvars < nnewlinvars);
            newlinvars  [inewlinvars] = var;
            newlincoefs[inewlinvars] = term->lincoefs;
            ++inewlinvars;
            SCIP_CALL( SCIPcaptureVar(scip, var) );
         }
      }
   assert(inewlinvars == nnewlinvars);
   consdata->nlinvars = nnewlinvars;
   consdata->linvars = newlinvars;
   consdata->lincoefs = newlincoefs;
   if( !SCIPisInfinity(scip, -consdata->lhs) )
      consdata->lhs -= constant;
   if( !SCIPisInfinity(scip, consdata->rhs) )
      consdata->rhs -= constant;
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->linrange, consdata->nlinvars) );

   SCIPdebugMessage("modified constraint %s to: ", SCIPconsGetName(cons));
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
#endif

   return SCIP_OKAY;
}
#endif

/** Reformulates products of binary variables as AND constraint.
 *  For a product y*x, with x and y binary variables, the product is replaced by a new auxiliary variable z and the constraint z = {x and y} is added.
 */
static
SCIP_RETCODE presolveTryAddAND(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_HASHMAP*         terms,              /**< constraint function in form of PRESOLVEQUADTERM's */
   int*                  nconsadded          /**< buffer where to add the number of AND constraints added */
   )
{
   SCIP_VAR*          x;
   SCIP_VAR*          y;
   PRESOLVEQUADTERM*  xterm;
   PRESOLVEQUADTERM*  yterm;
   PRESOLVEBILINITEM* bitem;
   char               name[255];
   SCIP_HASHMAPLIST*  list;
   SCIP_HASHMAPLIST*  blist;
   SCIP_VAR*          auxvar;
   SCIP_CONS*         andcons;
   SCIP_VAR*          vars[2];
   int                i;
   int                j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(terms != NULL);
   assert(nconsadded != NULL);
   
   *nconsadded = 0;
   
   for( i = 0; i < SCIPhashmapGetNLists(terms); ++i )
      for( list = SCIPhashmapGetList(terms, i); list != NULL; list = SCIPhashmapListGetNext(list) )
      {
         x = (SCIP_VAR*) SCIPhashmapListGetOrigin(list);
         if( SCIPvarGetType(x) != SCIP_VARTYPE_BINARY )
            continue;

         xterm = (PRESOLVEQUADTERM*) SCIPhashmapListGetImage(list);
         if( xterm->bilin == NULL )
            continue;

         for( j = 0; j < SCIPhashmapGetNLists(xterm->bilin); ++j )
         {
            blist = SCIPhashmapGetList(xterm->bilin, j);
            while( blist != NULL )
            {
               y = (SCIP_VAR*) SCIPhashmapListGetOrigin(blist);
               if( SCIPvarGetType(y) != SCIP_VARTYPE_BINARY )
               {
                  blist = SCIPhashmapListGetNext(blist);
                  continue;
               }

               (void)SCIPsnprintf(name, 255, "prod%s*%s", SCIPvarGetName(x), SCIPvarGetName(y));
               SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, 0., 1., 0., SCIP_VARTYPE_BINARY, TRUE, TRUE, NULL, NULL, NULL, NULL) );
               SCIP_CALL( SCIPaddVar(scip, auxvar) );

               vars[0] = x;
               vars[1] = y;
               (void)SCIPsnprintf(name, 255, "%sAND%s", SCIPvarGetName(x), SCIPvarGetName(y));
               SCIP_CALL( SCIPcreateConsAnd(scip, &andcons, name, auxvar, 2, vars,
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                  SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                  SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),   SCIPconsIsStickingAtNode(cons)) );
               SCIP_CALL( SCIPaddCons(scip, andcons) );
               SCIPdebugMessage("added AND constraint: ");
#ifdef SCIP_DEBUG
               SCIPprintCons(scip, andcons, NULL);
#endif
               SCIP_CALL( SCIPreleaseCons(scip, &andcons) );

               bitem = (PRESOLVEBILINITEM*) SCIPhashmapListGetImage(blist);
               SCIP_CALL( presolveQuadTermAdd(scip, terms, auxvar, bitem->coef, 0.0, NULL, 0.0) );
               SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

               /* delete BilinItem and entry from xterm->bilin */
               SCIPfreeBlockMemory(scip, &bitem);
               SCIP_CALL( SCIPhashmapRemove(xterm->bilin, y) );
               blist = SCIPhashmapGetList(xterm->bilin, j);

               /* find and delete BilinItem in yterm and entry from yterm->bilin */
               yterm = (PRESOLVEQUADTERM*) SCIPhashmapGetImage(terms, y);
               assert(yterm        != NULL);
               assert(yterm->bilin != NULL);
               bitem = (PRESOLVEBILINITEM*) SCIPhashmapGetImage(yterm->bilin, x);
               assert(bitem != NULL);
               SCIPfreeBlockMemory(scip, &bitem);
               SCIP_CALL( SCIPhashmapRemove(yterm->bilin, x) );

               ++*nconsadded;
            }
         }
      }

   if( *nconsadded != 0 )
   {
      /* clean empty bilin's and empty quad terms */
      for( i = 0; i < SCIPhashmapGetNLists(terms); ++i )
      {
         list = SCIPhashmapGetList(terms, i);
         while( list != NULL )
         {
            yterm = (PRESOLVEQUADTERM*) SCIPhashmapListGetImage(list);
            if( yterm->bilin != NULL && SCIPhashmapIsEmpty(yterm->bilin) )
               SCIPhashmapFree(&yterm->bilin);

            if( yterm->bilin == NULL && yterm->lincoefs == 0.0 && yterm->sqrcoefs == 0.0 )
            {
               SCIPfreeBlockMemory(scip, &yterm);
               SCIP_CALL( SCIPhashmapRemove(terms, (SCIP_VAR*)SCIPhashmapListGetOrigin(list)) );
               list = SCIPhashmapGetList(terms, i); /* restart from beginning of list */
            }
            else
               list = SCIPhashmapListGetNext(list);
         }
      }
   }

   return SCIP_OKAY;
}

/** Reformulates products of binary times bounded continuous variables as system of linear inequalities (plus auxiliary variable).
 * 
 *  For a product x*y, with y a binary variable and x a continous variable with finite bounds,
 *  an auxiliary variable z and the inequalities \f$ x^L * y \leq z \leq x^U * y \f$ and \f$ x - (1-y)*x^U \leq z \leq x - (1-y)*x^L \f$ are added.
 * 
 *  If x is a linear term consisting of more than one variable, it is split up in groups of linear terms of length at most maxnrvar.
 *  For each product of linear term of length at most maxnrvar with y, an auxiliary z and linear inequalities are added.
 * 
 *  If y is a binary variable, the AND constraint \f$ z = x \wedge y \f$ is added instead of linear constraints.
 */
static
SCIP_RETCODE presolveTryAddLinearReform(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_HASHMAP*         terms,              /**< constraint function in form of PRESOLVEQUADTERM's */
   int*                  nvarsadded,         /**< buffer where to store the number of auxiliary variables added */
   int                   maxnrvar,           /**< maximal number of variables in linear term to consider when replacing by one auxiliary variable */
   int                   empathy4and         /**< empathy for using AND constraint handler */
   )
{  /*lint --e{666} */
   SCIP_VAR**         xvars  = NULL;
   SCIP_Real*         xcoef = NULL;
   SCIP_INTERVAL      xbnds;
   SCIP_INTERVAL      tmp;
   int                nxvars;
   SCIP_VAR*          y;
   SCIP_VAR*          bvar;
   SCIP_HASHMAPLIST*  list;
   SCIP_HASHMAPLIST*  blist;
   PRESOLVEQUADTERM*  yterm;
   PRESOLVEBILINITEM* bitem;
   PRESOLVEQUADTERM*  bterm;
   char               name[255];
   int                nbilinterms;
   SCIP_VAR*          auxvar;
   SCIP_CONS*         auxcons;
   SCIP_Bool          maxnrvarfull; /* indicates whether we stopped collecting xvars because the maxnrvar limit was reached */
   int                i;
   int                j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(terms != NULL);
   assert(nvarsadded != NULL);
   
   *nvarsadded = 0;
   maxnrvarfull = FALSE;

   if( maxnrvar == 0 )
      return SCIP_OKAY;
   
   for( i = 0; i < SCIPhashmapGetNLists(terms); i = maxnrvarfull ? i : i+1 )
      for( list = SCIPhashmapGetList(terms, i); list != NULL; list = maxnrvarfull ? list : SCIPhashmapListGetNext(list) )
      {
         maxnrvarfull = FALSE;
         
         y = (SCIP_VAR*) SCIPhashmapListGetOrigin(list);
         if( SCIPvarGetType(y) != SCIP_VARTYPE_BINARY )
            continue;
         
         yterm = (PRESOLVEQUADTERM*) SCIPhashmapListGetImage(list);
         if( yterm->bilin == NULL )
            continue;
         
         nbilinterms = SCIPhashmapGetNEntries(yterm->bilin);
         if( nbilinterms == 0 )
            continue;
         
         nxvars = 0;
         SCIPintervalSet(&xbnds, 0.0);
         SCIP_CALL( SCIPreallocBufferArray(scip, &xvars,  MIN(maxnrvar, nbilinterms)+2) ); /* add 2 for later use when creating linear constraints */
         SCIP_CALL( SCIPreallocBufferArray(scip, &xcoef, MIN(maxnrvar, nbilinterms)+2) );
         
         /* setup list of variables x_i with coefficients a_i that are multiplied with binary y: y*(sum_i a_i*x_i)
          * and compute range of sum_i a_i*x_i
          */
         for( j = 0; j < SCIPhashmapGetNLists(yterm->bilin) && !maxnrvarfull; ++j )
            for( blist = SCIPhashmapGetList(yterm->bilin, j); blist != NULL; blist = SCIPhashmapListGetNext(blist) )
            {
               if( nxvars >= maxnrvar )
               {
                  maxnrvarfull = TRUE;
                  break;
               }
               
               bvar  = (SCIP_VAR*)          SCIPhashmapListGetOrigin(blist);
               if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(bvar)) || SCIPisInfinity(scip, SCIPvarGetUbGlobal(bvar)) )
                  continue;
               
               bitem = (PRESOLVEBILINITEM*) SCIPhashmapListGetImage(blist);
               xvars [nxvars] = bvar;
               xcoef[nxvars] = bitem->coef;
               
               SCIPintervalSetBounds(&tmp, MIN(SCIPvarGetLbGlobal(bvar), SCIPvarGetUbGlobal(bvar)), MAX(SCIPvarGetLbGlobal(bvar), SCIPvarGetUbGlobal(bvar)));
               SCIPintervalMulScalar(SCIPinfinity(scip), &tmp, tmp, bitem->coef);
               SCIPintervalAdd(SCIPinfinity(scip), &xbnds, xbnds, tmp);

               ++nxvars;
               
               /* free PRESOLVEBILINITEM structure in bilin map at y; entry from bilin map is removed later */
               SCIPfreeBlockMemory(scip, &bitem);

               /* free PRESOLVEBILINITEM structure in bilin map at bvar and delete entry from this bilin map */
               bterm = (PRESOLVEQUADTERM*)SCIPhashmapGetImage(terms, bvar);
               assert(bterm != NULL);
               bitem = (PRESOLVEBILINITEM*) SCIPhashmapGetImage(bterm->bilin, y);
               assert(bitem != NULL);
               SCIPfreeBlockMemory(scip, &bitem);
               SCIP_CALL( SCIPhashmapRemove(bterm->bilin, y) );
            }
         
         if( nxvars == 0 ) /* all x_i seem to be unbounded */
            continue;
         
         /* remove entries from bilin map at y */
         for( j = 0; j < nxvars; ++j )
         {
            SCIP_CALL( SCIPhashmapRemove(yterm->bilin, xvars[j]) );
         }
         
         assert(!SCIPisInfinity(scip, -SCIPintervalGetInf(xbnds)));
         assert(!SCIPisInfinity(scip, SCIPintervalGetSup(xbnds)));
         
         if( nxvars == 1 && empathy4and >= 1 && SCIPvarGetType(xvars[0]) == SCIP_VARTYPE_BINARY )
         { /* product of two binary variables, replace by auxvar and AND constraint */
            /* add auxiliary variable z */
            (void)SCIPsnprintf(name, 255, "prod%s*%s", SCIPvarGetName(y), SCIPvarGetName(xvars[0]));
            SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, 0., 1., 0., SCIP_VARTYPE_BINARY, TRUE, TRUE, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, auxvar) );
            
            xvars[1] = y;
            (void)SCIPsnprintf(name, 255, "%sAND%s", SCIPvarGetName(y), SCIPvarGetName(xvars[0]));
            SCIP_CALL( SCIPcreateConsAnd(scip, &auxcons, name, auxvar, 2, xvars,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
               SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),   SCIPconsIsStickingAtNode(cons)) );
            SCIP_CALL( SCIPaddCons(scip, auxcons) );
            SCIPdebugMessage("added AND constraint: ");
#ifdef SCIP_DEBUG
            SCIPprintCons(scip, auxcons, NULL);
#endif
            SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
            
            SCIP_CALL( presolveQuadTermAdd(scip, terms, auxvar, xcoef[0], 0.0, NULL, 0.0) );

            SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );
         }
         else
         { /* product of binary avariable with more than one binary or with continuous variables or with binary and user did not like AND -> replace by auxvar and linear constraints */
            /* add auxiliary variable z */
            if( nxvars == 1 )
               (void)SCIPsnprintf(name, 255, "prod%s*%s", SCIPvarGetName(y), SCIPvarGetName(xvars[0]));
            else
               (void)SCIPsnprintf(name, 255, "prod%s*%s*more", SCIPvarGetName(y), SCIPvarGetName(xvars[0]));
            SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, MIN(0., SCIPintervalGetInf(xbnds)), MAX(0., SCIPintervalGetSup(xbnds)), 0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, auxvar) );

            if( !SCIPisZero(scip, SCIPintervalGetInf(xbnds)) )
            { /* add 0 <= z - xbnds.inf * y constraint (as varbound constraint) */
               (void)SCIPsnprintf(name, 255, "linreform%s_1", SCIPvarGetName(y));
               SCIP_CALL( SCIPcreateConsVarbound(scip, &auxcons, name, auxvar, y, -SCIPintervalGetInf(xbnds), 0.0, SCIPinfinity(scip),
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                  SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                  SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),   SCIPconsIsStickingAtNode(cons)) );
               SCIP_CALL( SCIPaddCons(scip, auxcons) );
               SCIPdebugMessage("added varbound constraint: ");
#ifdef SCIP_DEBUG
               SCIPprintCons(scip, auxcons, NULL);
#endif
               SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
            }
            if( !SCIPisZero(scip, SCIPintervalGetSup(xbnds)) )
            { /* add z - xbnds.sup * y <= 0 constraint (as varbound constraint) */
               (void)SCIPsnprintf(name, 255, "linreform%s_2", SCIPvarGetName(y));
               SCIP_CALL( SCIPcreateConsVarbound(scip, &auxcons, name, auxvar, y, -SCIPintervalGetSup(xbnds), -SCIPinfinity(scip), 0.0,
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                  SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                  SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),   SCIPconsIsStickingAtNode(cons)) );
               SCIP_CALL( SCIPaddCons(scip, auxcons) );
               SCIPdebugMessage("added varbound constraint: ");
#ifdef SCIP_DEBUG
               SCIPprintCons(scip, auxcons, NULL);
#endif
               SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
            }

            /* add xbnds.inf <= sum_i a_i*x_i + xbnds.inf * y - z constraint */
            xvars[nxvars] = y;
            xvars[nxvars+1] = auxvar;
            xcoef[nxvars] = SCIPintervalGetInf(xbnds);
            xcoef[nxvars+1] = -1;

            (void)SCIPsnprintf(name, 255, "linreform%s_3", SCIPvarGetName(y));
            SCIP_CALL( SCIPcreateConsLinear(scip, &auxcons, name, nxvars+2, xvars, xcoef, SCIPintervalGetInf(xbnds), SCIPinfinity(scip),
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
               SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
            SCIP_CALL( SCIPaddCons(scip, auxcons) );
            SCIPdebugMessage("added linear constraint: ");
#ifdef SCIP_DEBUG
            SCIPprintCons(scip, auxcons, NULL);
#endif
            SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );

            /* add sum_i a_i*x_i + xbnds.sup * y - z <= xbnds.sup constraint */
            xcoef[nxvars] = SCIPintervalGetSup(xbnds);

            (void)SCIPsnprintf(name, 255, "linreform%s_4", SCIPvarGetName(y));
            SCIP_CALL( SCIPcreateConsLinear(scip, &auxcons, name, nxvars+2, xvars, xcoef, -SCIPinfinity(scip), SCIPintervalGetSup(xbnds),
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
               SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
            SCIP_CALL( SCIPaddCons(scip, auxcons) );
            SCIPdebugMessage("added linear constraint: ");
#ifdef SCIP_DEBUG
            SCIPprintCons(scip, auxcons, NULL);
#endif
            SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );

            /* add z to this constraint */
            SCIP_CALL( presolveQuadTermAdd(scip, terms, auxvar, 1.0, 0.0, NULL, 0.0) );
            
            SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );
         }
         
         ++*nvarsadded;
      }
      
   if( *nvarsadded != 0 )
   {
      SCIPdebugMessage("added %d sets of constraints to reformulate product with binary variable in constraint %s\n", *nvarsadded, SCIPconsGetName(cons));

      /* clean empty bilin's and empty quad terms */
      for( i = 0; i < SCIPhashmapGetNLists(terms); ++i )
      {
         list = SCIPhashmapGetList(terms, i);
         while( list != NULL )
         {
            yterm = (PRESOLVEQUADTERM*) SCIPhashmapListGetImage(list);
            if( yterm->bilin != NULL && SCIPhashmapIsEmpty(yterm->bilin) )
               SCIPhashmapFree(&yterm->bilin);
            
            if( yterm->bilin == NULL && yterm->lincoefs == 0.0 && yterm->sqrcoefs == 0.0 )
            {
               SCIPfreeBlockMemory(scip, &yterm);
               SCIP_CALL( SCIPhashmapRemove(terms, (SCIP_VAR*)SCIPhashmapListGetOrigin(list)) );
               list = SCIPhashmapGetList(terms, i); /* restart from beginning of list */
            }
            else
               list = SCIPhashmapListGetNext(list);
         }
      }
   }

   SCIPfreeBufferArrayNull(scip, &xvars);
   SCIPfreeBufferArrayNull(scip, &xcoef);

   return SCIP_OKAY;
}

/** tries to automatically convert a quadratic constraint (or a part of it) into a more specific and more specialized constraint */
static
SCIP_RETCODE presolveUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler data structure */
   SCIP_CONS*            cons,               /**< source constraint to try to convert */
   SCIP_CONS**           upgdconslhs,        /**< pointer to store upgraded constraint for left  hand side, or NULL if not successful */
   SCIP_CONS**           upgdconsrhs         /**< pointer to store upgraded constraint for right hand side, or NULL if not successful */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* localupgdconslhs;
   SCIP_CONS* localupgdconsrhs;
   SCIP_VAR* var;
   SCIP_Real lincoef;
   SCIP_Real quadcoef;
   SCIP_Real lb;
   SCIP_Real ub;
   int nbinlin;
   int nbinquad;
   int nintlin;
   int nintquad;
   int nimpllin;
   int nimplquad;
   int ncontlin;
   int ncontquad;
   SCIP_Bool integral;
   int i;

   assert(upgdconslhs != NULL);
   assert(upgdconsrhs != NULL);

   *upgdconslhs = NULL;
   *upgdconsrhs = NULL;

   /* we cannot upgrade a modifiable linear constraint, since we don't know what additional coefficients to expect */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   /* get constraint handler data and constraint data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   consdata     = SCIPconsGetData(cons);
   assert(conshdlrdata != NULL);
   assert(consdata     != NULL);

   /*
    * calculate some statistics on quadratic constraint
    */

   nbinlin   = 0;
   nbinquad  = 0;
   nintlin   = 0;
   nintquad  = 0;
   nimpllin  = 0;
   nimplquad = 0;
   ncontlin  = 0;
   ncontquad = 0;
   integral  = TRUE;
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      var = consdata->linvars[i];
      lincoef = consdata->lincoefs[i];
      lb  = SCIPvarGetLbLocal(var);
      ub  = SCIPvarGetUbLocal(var);
      assert(!SCIPisZero(scip, lincoef));

      switch( SCIPvarGetType(var) )
      {
      case SCIP_VARTYPE_BINARY:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef);
         nbinlin++;
         break;
      case SCIP_VARTYPE_INTEGER:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef);
         nintlin++;
         break;
      case SCIP_VARTYPE_IMPLINT:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef);
         nimpllin++;
         break;
      case SCIP_VARTYPE_CONTINUOUS:
         integral = integral && SCIPisEQ(scip, lb, ub) && SCIPisIntegral(scip, lincoef * lb);
         ncontlin++;
         break;
      default:
         SCIPerrorMessage("unknown variable type\n");
         return SCIP_INVALIDDATA;
      }
   }

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      var = consdata->quadvars[i];
      lincoef  = consdata->quadlincoefs[i];
      quadcoef = consdata->quadlincoefs[i];
      lb  = SCIPvarGetLbLocal(var);
      ub  = SCIPvarGetUbLocal(var);

      switch( SCIPvarGetType(var) )
      {
      case SCIP_VARTYPE_BINARY:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef) && SCIPisIntegral(scip, quadcoef);
         nbinquad++;
         break;
      case SCIP_VARTYPE_INTEGER:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef) && SCIPisIntegral(scip, quadcoef);
         nintquad++;
         break;
      case SCIP_VARTYPE_IMPLINT:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef) && SCIPisIntegral(scip, quadcoef);
         nimplquad++;
         break;
      case SCIP_VARTYPE_CONTINUOUS:
         integral = integral && SCIPisEQ(scip, lb, ub) && SCIPisIntegral(scip, lincoef * lb + quadcoef * lb * lb);
         ncontquad++;
         break;
      default:
         SCIPerrorMessage("unknown variable type\n");
         return SCIP_INVALIDDATA;
      }
   }
   
   if( integral )
   {
      for( i = 0; i < consdata->nbilinterms && integral; ++i )
      {
         if( SCIPvarGetType(consdata->bilinvars1[i]) < SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(consdata->bilinvars2[i]) < SCIP_VARTYPE_CONTINUOUS )
            integral = integral && SCIPisIntegral(scip, consdata->bilincoefs[i]);
         else
            integral = FALSE;
      }
   }

   /*
    * call the upgrading methods
    */

   SCIPdebugMessage("upgrading quadratic constraint <%s> (%d upgrade methods):\n",
      SCIPconsGetName(cons), conshdlrdata->nquadconsupgrades);
   SCIPdebugMessage(" binlin=%d binquad=%d intlin=%d intquad=%d impllin=%d implquad=%d contlin=%d contquad=%d integral=%u\n",
      nbinlin, nbinquad, nintlin, nintquad, nimpllin, nimplquad, ncontlin, ncontquad, integral);

   /* try all upgrading methods in priority order in case the upgrading step is enable  */
   for( i = 0; i < conshdlrdata->nquadconsupgrades; ++i )
   {
      if( conshdlrdata->quadconsupgrades[i]->active )
      {
         localupgdconslhs = NULL;
         localupgdconsrhs = NULL;
         
         SCIP_CALL( conshdlrdata->quadconsupgrades[i]->quadconsupgd(scip, cons,
            nbinlin, nbinquad, nintlin, nintquad, nimpllin, nimplquad, ncontlin, ncontquad, integral,
            &localupgdconslhs, &localupgdconsrhs) );

         if( localupgdconslhs != NULL )
         { /* got upgrade for left hand side */
            assert(!SCIPisInfinity(scip, -consdata->lhs));
            SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
            SCIPdebugMessage(" -> upgraded to constraint type <%s>\n", SCIPconshdlrGetName(SCIPconsGetHdlr(localupgdconslhs)));
            SCIPdebug( SCIP_CALL( SCIPprintCons(scip, localupgdconslhs, NULL) ) );
            
            *upgdconslhs  = localupgdconslhs;
            consdata->lhs = -SCIPinfinity(scip);
         }
         
         if( localupgdconsrhs != NULL )
         { /* got upgrade for right hand side, but maybe same constraint as lhs */
            assert(!SCIPisInfinity(scip,  consdata->rhs));
            if (localupgdconslhs != localupgdconsrhs)
            { /* upgrade for rhs is not the same constraint as for lhs */
               SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
               SCIPdebugMessage(" -> upgraded to constraint type <%s>\n", SCIPconshdlrGetName(SCIPconsGetHdlr(localupgdconsrhs)));
               SCIPdebug( SCIP_CALL( SCIPprintCons(scip, localupgdconsrhs, NULL) ) );
            }
            
            *upgdconsrhs  = localupgdconsrhs;
            consdata->rhs = SCIPinfinity(scip);
         }
         
         /* both sides of constraint have been upgraded, so we can stop here */
         if( SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, consdata->rhs) )
            break;
      }
   }

   return SCIP_OKAY;
}

/** checks a quadratic constraint for convexity and/or concavity */
static
SCIP_RETCODE checkCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< quadratic constraint handler */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;
   double*        matrix;
   SCIP_HASHMAP*  var2index;
   int            i;
   int            n;
   int            nn;
#ifdef WITH_LAPACK
   int            row;
   int            col;
   double*        alleigval;
#endif

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   n = consdata->nquadvars;

   if( n == 1 )
   {
      assert(consdata->nbilinterms == 0);
      consdata->isconvex = !SCIPisNegative(scip, consdata->quadsqrcoefs[0]);
      consdata->isconcave = !SCIPisPositive(scip, consdata->quadsqrcoefs[0]);
      return SCIP_OKAY;
   }

   if( n == 0 )
   {
      consdata->isconvex = TRUE;
      consdata->isconcave = TRUE;
      return SCIP_OKAY;
   }

   if( consdata->nbilinterms == 0 )
   {
      consdata->isconvex = TRUE;
      consdata->isconcave = TRUE;
      for( i = 0; i < n; ++i )
      {
         consdata->isconvex  = consdata->isconvex  && !SCIPisNegative(scip, consdata->quadsqrcoefs[i]);
         consdata->isconcave = consdata->isconcave && !SCIPisPositive(scip, consdata->quadsqrcoefs[i]);
      }
      return SCIP_OKAY;
   }

   if( n == 2 )
   { /* compute eigenvalues by hand */
      assert(consdata->nbilinterms == 1);
      consdata->isconvex =
         consdata->quadsqrcoefs[0] >= 0 &&
         consdata->quadsqrcoefs[1] >= 0 &&
         4 * consdata->quadsqrcoefs[0] * consdata->quadsqrcoefs[1] >= consdata->bilincoefs[0] * consdata->bilincoefs[0];
      consdata->isconcave = 
         consdata->quadsqrcoefs[0] <= 0 &&
         consdata->quadsqrcoefs[1] <= 0 &&
         4 * consdata->quadsqrcoefs[0] * consdata->quadsqrcoefs[1] >= consdata->bilincoefs[0] * consdata->bilincoefs[0];
      return SCIP_OKAY;
   }

   /* lower triangular of quadratic term matrix, scaled by box diameter */
   nn = n * n;
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, nn) );
   BMSclearMemoryArray(matrix, nn);

   consdata->isconvex  = TRUE;
   consdata->isconcave = TRUE;

   SCIP_CALL( SCIPhashmapCreate(&var2index, SCIPblkmem(scip), n) );
   for( i = 0; i < n; ++i )
   {
      if( consdata->nadjbilin[i] )
      {
         SCIP_CALL( SCIPhashmapInsert(var2index, consdata->quadvars[i], (void*)(size_t)i) );
         matrix[i*n + i] = consdata->quadsqrcoefs[i];
      }
      /* nonzero elements on diagonal tell a lot about convexity/concavity */
      if( SCIPisNegative(scip, consdata->quadsqrcoefs[i]) )
         consdata->isconvex  = FALSE;
      if( SCIPisPositive(scip, consdata->quadsqrcoefs[i]) )
         consdata->isconcave = FALSE;
   }

   if( !consdata->isconvex && !consdata->isconcave )
   {
      SCIPfreeBufferArray(scip, &matrix);
      SCIPhashmapFree(&var2index);
      return SCIP_OKAY;
   }

#ifdef WITH_LAPACK
   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      assert(SCIPhashmapExists(var2index, consdata->bilinvars1[i]));
      assert(SCIPhashmapExists(var2index, consdata->bilinvars2[i]));
      row = (int)(size_t)SCIPhashmapGetImage(var2index, consdata->bilinvars1[i]);
      col = (int)(size_t)SCIPhashmapGetImage(var2index, consdata->bilinvars2[i]);
      if( row < col )
         matrix[row * n + col] = consdata->bilincoefs[i]/2;
      else
         matrix[col * n + row] = consdata->bilincoefs[i]/2;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &alleigval, n) );
   /* TODO can we compute only min and max eigval?
      TODO can we estimate the numerical error? */
   if( LapackDsyev(scip, FALSE, n, matrix, alleigval) != SCIP_OKAY )
   {
      SCIPwarningMessage("Failed to compute eigenvalues of quadratic coefficient matrix of constraint %s. Assuming matrix is indefinite.\n", SCIPconsGetName(cons));
      consdata->isconvex = FALSE;
      consdata->isconcave = FALSE;
   }
   else
   {
      consdata->isconvex &= !SCIPisNegative(scip, alleigval[0]);
      consdata->isconcave &= !SCIPisPositive(scip, alleigval[n-1]);
   }
   SCIPfreeBufferArray(scip, &alleigval);
   
#else
   consdata->isconvex = FALSE;
   consdata->isconcave = FALSE;
#endif
   
   SCIPhashmapFree(&var2index);
   SCIPfreeBufferArray(scip, &matrix);

   return SCIP_OKAY;
}

/** locks a linear variable in a constraint */
static
SCIP_RETCODE lockLinearVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where to lock a variable */
   SCIP_VAR*             var,                /**< variable to lock */
   SCIP_Real             coef                /**< coefficient of variable in constraint */
   )
{
   SCIP_CONSDATA* consdata;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != 0.0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( coef > 0.0 )
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip,  consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons, !SCIPisInfinity(scip,  consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** unlocks a linear variable in a constraint */
static
SCIP_RETCODE unlockLinearVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where to unlock a variable */
   SCIP_VAR*             var,                /**< variable to unlock */
   SCIP_Real             coef                /**< coefficient of variable in constraint */
   )
{
   SCIP_CONSDATA* consdata;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != 0.0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( coef > 0.0 )
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip,  consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, !SCIPisInfinity(scip,  consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** sets bounds for variables in not evidently convex terms to some predefined value */
static
SCIP_RETCODE boundUnboundedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real             bound,              /**< value to use for bound */
   int*                  nchgbnds            /**< buffer where to add the number of bound changes, or NULL */
   )
{
   SCIP_Bool      infeasible;
   SCIP_CONSDATA* consdata;
   int            i;
   
   assert(scip != NULL);
   assert(cons != NULL);

   if( SCIPisInfinity(scip, bound) )
      return SCIP_OKAY;

   consdata =  SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( consdata->nadjbilin[i] == 0 &&
          (SCIPisInfinity(scip,  consdata->rhs) || consdata->quadsqrcoefs[i] > 0) &&
          (SCIPisInfinity(scip, -consdata->lhs) || consdata->quadsqrcoefs[i] < 0) )
         continue; /* skip evidently convex terms */

      if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvars[i])) )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "set lower bound of %s to %g\n", SCIPvarGetName(consdata->quadvars[i]), -bound);
         SCIP_CALL( SCIPtightenVarLb(scip, consdata->quadvars[i], -bound, FALSE, &infeasible, NULL) );
         assert(!infeasible);
         if( nchgbnds != NULL )
            ++*nchgbnds;
      }

      if( SCIPisInfinity(scip,  SCIPvarGetUbLocal(consdata->quadvars[i])) )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "set upper bound of %s to %g\n", SCIPvarGetName(consdata->quadvars[i]),  bound);
         SCIP_CALL( SCIPtightenVarUb(scip, consdata->quadvars[i],  bound, FALSE, &infeasible, NULL) );
         assert(!infeasible);
         if( nchgbnds != NULL ) 
            ++*nchgbnds;
      }
   }

   return SCIP_OKAY;
}

/** gets euclidean norm of gradient of quadratic function */
static
SCIP_Real getGradientNorm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol                 /**< solution or NULL if LP solution should be used */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      norm = 0.0;
   SCIP_Real      g;
   int            i, j, k;
   SCIP_VAR*      var;
   
   assert(scip != NULL);
   assert(cons != NULL);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   /* TODO allow also other norms than euclidean, maybe read separating/efficacynorm */
   
   for( i = 0; i < consdata->nlinvars; ++i )
      norm += consdata->lincoefs[i] * consdata->lincoefs[i];
   
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      var = consdata->quadvars[i];
      assert(!SCIPisInfinity(scip,  SCIPgetSolVal(scip, sol, var)));
      assert(!SCIPisInfinity(scip, -SCIPgetSolVal(scip, sol, var)));
      g  =     consdata->quadlincoefs[i];
      g += 2 * consdata->quadsqrcoefs[i] * SCIPgetSolVal(scip, sol, var);
      for( j = 0; j < consdata->nadjbilin[i]; ++j )
      {
         k = consdata->adjbilin[i][j];
         if( consdata->bilinvars1[k] == var )
            g += consdata->bilincoefs[k] * SCIPgetSolVal(scip, sol, consdata->bilinvars2[k]);
         else
            g += consdata->bilincoefs[k] * SCIPgetSolVal(scip, sol, consdata->bilinvars1[k]);
      }
      norm += g*g;
   }
   
   return sqrt(norm);
}

/** computes violation of a constraint */
static
SCIP_RETCODE computeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   SCIP_Bool             doscaling           /**< should we scale the violation by the gradient of the quadratic function ? */ 
   )
{  /*lint --e{666}*/
   SCIP_CONSDATA* consdata;
   SCIP_Real val;
   SCIP_Real varval;
   int i;
   int j;
   
   assert(scip != NULL);
   assert(cons != NULL);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   val = 0.0;
   varval = 0.0;
   
   /* TODO take better care of variables at +/- infinity: e.g., run instance waste in debug mode with a short timelimit (30s) */
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      if( SCIPisInfinity(scip, ABS(SCIPgetSolVal(scip, sol, consdata->linvars[i]))) )
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            consdata->lhsviol = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            consdata->rhsviol = SCIPinfinity(scip);
         return SCIP_OKAY;
      }
      val += consdata->lincoefs[i] * SCIPgetSolVal(scip, sol, consdata->linvars[i]);
   }

   for( j = 0; j < consdata->nquadvars; ++j )
   {
      varval = SCIPgetSolVal(scip, sol, consdata->quadvars[j]);
      if( SCIPisInfinity(scip, ABS(varval)) )
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            consdata->lhsviol = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            consdata->rhsviol = SCIPinfinity(scip);
         return SCIP_OKAY;
      }
      val   += (consdata->quadlincoefs[j] + consdata->quadsqrcoefs[j] * varval) * varval;
   }
   
   for( j = 0; j < consdata->nbilinterms; ++j )
      val += consdata->bilincoefs[j] * SCIPgetSolVal(scip, sol, consdata->bilinvars1[j]) * SCIPgetSolVal(scip, sol, consdata->bilinvars2[j]);

   if( val < consdata->lhs && !SCIPisInfinity(scip, -consdata->lhs) )
      consdata->lhsviol = consdata->lhs - val;
   else
      consdata->lhsviol = 0.0;
   
   if( val > consdata->rhs && !SCIPisInfinity(scip,  consdata->rhs) )
      consdata->rhsviol = val - consdata->rhs;
   else
      consdata->rhsviol = 0.0;
   
   if( doscaling && (consdata->lhsviol || consdata->rhsviol) )
   {
      SCIP_Real norm = getGradientNorm(scip, cons, sol);
      if( norm > 1.0 )
      { /* TODO scale only if > 1., or should it be larger SCIPsumepsilon? */
         consdata->lhsviol /= norm;
         consdata->rhsviol /= norm;
      }
   }
   
   return SCIP_OKAY;
}

/** computes violation of a set of constraints */
static
SCIP_RETCODE computeViolations(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   SCIP_Bool             doscaling,          /**< are we scaling when computing violation ? */
   SCIP_CONS**           maxviolcon          /**< buffer to store constraint with largest violation, or NULL if solution is feasible */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      viol;
   SCIP_Real      maxviol;
   int            c;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(maxviolcon != NULL);
   
   *maxviolcon = NULL;

   maxviol = 0.0;
   
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);
      
      SCIP_CALL( computeViolation(scip, conss[c], sol, doscaling) );

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      
      viol = MAX(consdata->lhsviol, consdata->rhsviol);
      if( viol > maxviol && SCIPisFeasPositive(scip, viol) )
      {
         maxviol = viol;
         *maxviolcon = conss[c];
      }
   }
   
   return SCIP_OKAY;
}

/** Adds range of quadratic term w.r.t. local bounds to given interval.
 */
static
SCIP_RETCODE addQuadRange(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real             intervalinfty,      /**< value of infinity to use for interval operations */
   SCIP_INTERVAL*        resultant,          /**< interval where to add to */
   SCIP_VAR*             except              /**< a variable to skip in evaluation, NULL if nothing should be skipped */
   )
{  /*lint --e{666}*/
   SCIP_CONSDATA* consdata;
   SCIP_INTERVAL  lincoefs; /* linear interval coefficient in quadratic form */
   SCIP_INTERVAL  xrng;
   SCIP_INTERVAL  tmp;
   int            i;
   int            j; 
   int            k;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(resultant != NULL);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( consdata->quadvars[i] == except )
         continue;
      
      SCIPintervalSetBounds(&xrng, 
         -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->quadvars[i]), SCIPvarGetUbLocal(consdata->quadvars[i]))),
          infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->quadvars[i]), SCIPvarGetUbLocal(consdata->quadvars[i]))));
      SCIPintervalSet(&lincoefs, consdata->quadlincoefs[i]);

      for( j = 0; j < consdata->nadjbilin[i]; ++j )
      {
         k = consdata->adjbilin[i][j];
         if( consdata->bilinvars1[k] != consdata->quadvars[i] )
            continue; /* handle this term later */
         if( consdata->bilinvars2[k] == except )
            continue; /* variable is skipped */
         
         SCIPintervalSetBounds(&tmp, 
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvars2[k]), SCIPvarGetUbLocal(consdata->bilinvars2[k]))),
             infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvars2[k]), SCIPvarGetUbLocal(consdata->bilinvars2[k]))));
         SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilincoefs[k]);
         SCIPintervalAdd(intervalinfty, &lincoefs, lincoefs, tmp);
      }
      
      SCIPintervalQuad(intervalinfty, &tmp, consdata->quadsqrcoefs[i], lincoefs, xrng);
      assert(SCIPintervalGetSup(tmp) > -intervalinfty);
      assert(SCIPintervalGetInf(tmp) <  intervalinfty);
      SCIPintervalAdd(intervalinfty, resultant, *resultant, tmp);
   }
   
   return SCIP_OKAY;
}

/** checks by interval analysis whether a violated constraint is infeasible
 * 
 *  If lhsviol and rhsviol is below feasibility tolerance, the check is skipped.
 * 
 *  If isfeasible is set to false, then constraint is infeasible w.r.t. current local bounds.
 * 
 *  If isfeasible is set to true, then this gives no information.
 */
static
SCIP_RETCODE isIntervalFeasible(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_Bool*            isfeasible          /**< buffer to store the result */
   )
{  /*lint --e{666}*/
   SCIP_CONSDATA* consdata;
   SCIP_INTERVAL  val;
   SCIP_Real      intervalinfty;
   int            i;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(isfeasible != NULL);
   assert(SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVING && SCIPgetStage(scip) < SCIP_STAGE_SOLVED);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   *isfeasible = TRUE;
   intervalinfty = 1000 * SCIPinfinity(scip) * SCIPinfinity(scip);
   
   if( !SCIPisPositive(scip, consdata->lhsviol) && !SCIPisPositive(scip, consdata->rhsviol) )
      return SCIP_OKAY; /* obviously we have a feasible point */

   if( SCIPintervalIsEmpty(consdata->quadrange) )
   { /* need to update quadrange */
      SCIPintervalSet(&consdata->quadrange, 0.0);
      SCIP_CALL( addQuadRange(scip, cons, intervalinfty, &consdata->quadrange, NULL) );
   }

   val = consdata->quadrange;
   
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      if( SCIPintervalIsEmpty(consdata->linrange[i]) )
      { /* need to update linrange for var. i */
         SCIPintervalSetBounds(&consdata->linrange[i], 
            -infty2infty(SCIPinfinity(scip), intervalinfty, -SCIPvarGetLbLocal(consdata->linvars[i])), 
             infty2infty(SCIPinfinity(scip), intervalinfty,  SCIPvarGetUbLocal(consdata->linvars[i])));
         SCIPintervalMulScalar(intervalinfty, &consdata->linrange[i], consdata->linrange[i], consdata->lincoefs[i]);
      }
      SCIPintervalAdd(intervalinfty, &val, val, consdata->linrange[i]);
   }
   
   if( SCIPisFeasGT(scip, consdata->lhs, SCIPintervalGetSup(val)) || SCIPisFeasLT(scip, consdata->rhs, SCIPintervalGetInf(val)) )
   {
      SCIPdebugMessage("interval arithmetic found constraint %s infeasible: bounds = [%g, %g], interval = [%g, %g]\n", SCIPconsGetName(cons), consdata->lhs, consdata->rhs, SCIPintervalGetInf(val), SCIPintervalGetSup(val));
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *isfeasible = FALSE;
   }
   
   return SCIP_OKAY;
}

/** checks by interval analysis whether a set of constraints is infeasible
 *  If isfeasible is set to false, then one constraint is infeasible w.r.t. current local bounds.
 *  If isfeasible is set to true, then this gives no information.
 */
static
SCIP_RETCODE areIntervalFeasible(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints to check */
   int                   nconss,             /**< number of constraints to check */
   SCIP_CONS*            firstcons,          /**< constraint to check first, can be NULL */
   SCIP_Bool*            isfeasible          /**< buffer to store the result */
   )
{
   int c;
   
   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(isfeasible != NULL);

   if( firstcons != NULL )
   {
      SCIP_CALL( isIntervalFeasible(scip, firstcons, isfeasible) );
      if( !*isfeasible )
         return SCIP_OKAY;
   }

   for( c = 0; *isfeasible && c < nconss; ++c )
   {
      assert(conss != NULL);
      if( conss[c] == firstcons )
         continue;
      SCIP_CALL( isIntervalFeasible(scip, conss[c], isfeasible) );
   }

   return SCIP_OKAY;
}

/** generates a cut based on linearization for a carefully choosen reference point
 */
static
SCIP_RETCODE generateCutCareful(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   SCIP_BOUNDTYPE        violbound,          /**< for which bound a cut should be generated */
   SCIP_ROW**            row,                /**< storage for cut */
   SCIP_Real             maxrange            /**< maximal range allowed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool      isconvex;
   SCIP_Bool      isglobal;
   SCIP_Real      coef;
   SCIP_Real      rowcoef;
   SCIP_Real      bnd;
   SCIP_VAR*      x;
   SCIP_VAR*      y;
   SCIP_Real      xval;
   SCIP_Real      yval;
   int            j;
   SCIP_HASHMAP*  quadvaridx = NULL;
   SCIP_Real*     refpoint;
   SCIP_Real      threshold;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row  != NULL);
   assert(*row == NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   isconvex = (violbound == SCIP_BOUNDTYPE_LOWER) ? consdata->isconcave : consdata->isconvex;
   isglobal = SCIPconsIsGlobal(cons) && isconvex;
   
   assert(isconvex);
   
   /* create mapping quadvars to index in quadvars arrays, if have bilin terms */
   if( consdata->nbilinterms )
   {
      SCIP_CALL( SCIPhashmapCreate(&quadvaridx, SCIPblkmem(scip), consdata->nquadvars) );
      for( j = 0; j < consdata->nquadvars; ++j)
      {
         x = consdata->quadvars[j];
         SCIP_CALL( SCIPhashmapInsert(quadvaridx, x, (void*)(size_t)j) );
      }
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &refpoint, consdata->nquadvars) );
   
   for( threshold = 1E+6; threshold > 0.1; threshold /= 10 )
   {
      SCIPdebugMessage("try with threshold %g\n", threshold);
      assert(*row == NULL);
      
      SCIP_CALL( SCIPcreateEmptyRow(scip, row, "cut", -SCIPinfinity(scip), SCIPinfinity(scip), !isglobal /* locally */, FALSE /* modifiable */, TRUE /* removable */ ) );
      bnd = (violbound == SCIP_BOUNDTYPE_LOWER) ? consdata->lhs : consdata->rhs;
      assert(!SCIPisInfinity(scip, ABS(bnd)));

      /* add linear part */
      SCIP_CALL( SCIPaddVarsToRow(scip, *row, consdata->nlinvars, consdata->linvars, consdata->lincoefs) );

      /* create refpoint */
      for( j = 0; j < consdata->nquadvars; ++j)
      {
         x = consdata->quadvars[j];
         xval = SCIPgetSolVal(scip, sol, x);
         if( xval > threshold )
            xval = MAX( threshold, SCIPvarGetLbLocal(x));
         else if( xval < -threshold )
            xval = MIN(-threshold, SCIPvarGetUbLocal(x));

         refpoint[j] = xval;
      }

      for( j = 0; j < consdata->nquadvars; ++j )
      { /* linear term + linearization of square term */
         x = consdata->quadvars[j];
         rowcoef = consdata->quadlincoefs[j];

         if( consdata->quadsqrcoefs[j] )
         {
            xval = refpoint[j];

            rowcoef += 2*consdata->quadsqrcoefs[j]*xval;
            bnd     +=   consdata->quadsqrcoefs[j]*xval*xval;

            if( SCIPisInfinity(scip, ABS(rowcoef)) )
            {
               SCIPdebugMessage("skip linearization of square term in constraint %s because coeff. of var %s is at infinity (%g)\n", SCIPconsGetName(cons), SCIPvarGetName(x), rowcoef);
               SCIP_CALL( SCIPreleaseRow(scip, row) );
               break;
            }
         }

         SCIP_CALL( SCIPaddVarToRow(scip, *row, x, rowcoef) );
      }
      
      if( *row == NULL )
         continue;

      for( j = 0; j < consdata->nbilinterms; ++j )
      { /* linearization of bilinear terms */
         coef = consdata->bilincoefs[j];
         x    = consdata->bilinvars1[j];

         assert(SCIPhashmapExists(quadvaridx, x));
         assert((int)(size_t)SCIPhashmapGetImage(quadvaridx, x) < consdata->nquadvars);

         xval = refpoint[(int)(size_t)SCIPhashmapGetImage(quadvaridx, x)];

         if( SCIPisInfinity(scip, ABS(coef*xval)) )
         {
            SCIPdebugMessage("skip linearization of bilinear term in constraint %s because var %s gives coeff. at infinity (%g)\n", SCIPconsGetName(cons), SCIPvarGetName(x), coef*xval);
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            break;
         }

         y = consdata->bilinvars2[j];

         assert(SCIPhashmapExists(quadvaridx, y));
         assert((int)(size_t)SCIPhashmapGetImage(quadvaridx, y) < consdata->nquadvars);

         yval = refpoint[(int)(size_t)SCIPhashmapGetImage(quadvaridx, y)];

         if( SCIPisInfinity(scip, ABS(coef*yval)) )
         {
            SCIPdebugMessage("skip linearization of bilinear term in constraint %s because var %s gives coeff. at infinity (%g)\n", SCIPconsGetName(cons), SCIPvarGetName(y), coef*yval);
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            break;
         }

         SCIP_CALL( SCIPaddVarToRow(scip, *row, x, coef * yval) );
         SCIP_CALL( SCIPaddVarToRow(scip, *row, y, coef * xval) );
         bnd += coef * xval * yval;
      }
      
      if( *row == NULL )
         continue;

      SCIPdebugMessage(" -> generated careful cut rhs=%g, min=%f, max=%f range=%g\n",
          ABS(bnd),
          SCIPgetRowMinCoef(scip, *row), SCIPgetRowMaxCoef(scip, *row),
          SCIPgetRowMaxCoef(scip, *row)/SCIPgetRowMinCoef(scip, *row));

      if( SCIPisInfinity(scip, ABS(bnd)) )
      { /* seems to be a numerically bad cut */
         SCIPdebugMessage("skip cut for constraint %s because of very large left or right hand side: %g\n", SCIPconsGetName(cons), bnd);
         SCIP_CALL( SCIPreleaseRow(scip, row) );
         continue;
      }

      if( SCIPisGT(scip, SCIPgetRowMaxCoef(scip, *row)/SCIPgetRowMinCoef(scip, *row), maxrange) )
      { /* seems to be a numerically bad cut */
         SCIPdebugMessage("skip cut for constraint %s because of very large range: %g\n", SCIPconsGetName(cons), SCIPgetRowMaxCoef(scip, *row)/SCIPgetRowMinCoef(scip, *row));
         SCIP_CALL( SCIPreleaseRow(scip, row) );
         continue;
      }

      if( violbound == SCIP_BOUNDTYPE_LOWER )
      {
         SCIP_CALL( SCIPchgRowLhs(scip, *row, bnd) );
      }
      else
      {
         SCIP_CALL( SCIPchgRowRhs(scip, *row, bnd) );
      }
      break; /* have a cut that is good enough */
   }

   SCIPfreeBufferArray(scip, &refpoint);
   if( quadvaridx != NULL )
      SCIPhashmapFree(&quadvaridx);
   
   return SCIP_OKAY;
}

/** generates a cut based on linearization (if convex) or McCormick (if nonconvex)
 */
static
SCIP_RETCODE generateCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   SCIP_BOUNDTYPE        violbound,          /**< for which bound a cut should be generated */
   SCIP_ROW**            row,                /**< storage for cut */
   SCIP_Real             maxrange            /**< maximal range allowed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool      isconvex;
   SCIP_Bool      isglobal;
   SCIP_Real      coef;
   SCIP_Real      rowcoef;
   SCIP_Real      xcoef;
   SCIP_Real      ycoef;
   SCIP_Real      bnd;
   SCIP_Real      bnd_;
   SCIP_VAR*      x;
   SCIP_VAR*      y;
   SCIP_Real      xlb;
   SCIP_Real      xub;
   SCIP_Real      xval;
   SCIP_Real      ylb;
   SCIP_Real      yub;
   SCIP_Real      yval;
   int            j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row != NULL);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   isconvex = (violbound == SCIP_BOUNDTYPE_LOWER) ? consdata->isconcave : consdata->isconvex;
   isglobal = SCIPconsIsGlobal(cons) && isconvex;

   SCIP_CALL( SCIPcreateEmptyRow(scip, row, "cut", -SCIPinfinity(scip), SCIPinfinity(scip), !isglobal /* locally */, FALSE /* modifiable */, TRUE /* removable */ ) );
   bnd = (violbound == SCIP_BOUNDTYPE_LOWER) ? consdata->lhs : consdata->rhs;
   assert(!SCIPisInfinity(scip, ABS(bnd)));

   /* add linear part */
   SCIP_CALL( SCIPaddVarsToRow(scip, *row, consdata->nlinvars, consdata->linvars, consdata->lincoefs) );
   /* TODO should we buffer the coefficients of the quadratic variables before adding them to the row? */

   if( isconvex )
   {  /* do first-order taylor for each term */
      for( j = 0; j < consdata->nquadvars; ++j )
      { /* linear term + linearization of square term */
         x = consdata->quadvars[j];
         rowcoef = consdata->quadlincoefs[j];
         
         if( consdata->quadsqrcoefs[j] )
         {
            xval = SCIPgetSolVal(scip, sol, x);
            /* can happen when called from initlp */
            if( xval < SCIPvarGetLbLocal(x) )
               xval = SCIPvarGetLbLocal(x);
            else if( xval > SCIPvarGetUbLocal(x) )
               xval = SCIPvarGetUbLocal(x);
            if( SCIPisInfinity(scip, ABS(xval)) )
            {
               SCIPdebugMessage("skip linearization of square term in constraint %s because var %s is at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(x));
               SCIP_CALL( SCIPreleaseRow(scip, row) );
               return SCIP_OKAY;
            }

            if( consdata->nbilinterms || SCIPvarGetType(x) == SCIP_VARTYPE_CONTINUOUS || SCIPisIntegral(scip, xval) )
            {
               rowcoef += 2*consdata->quadsqrcoefs[j]*xval;
               bnd += consdata->quadsqrcoefs[j]*xval*xval;
            }
            else
            { /* if variable is discrete but fractional and there are no bilinear terms, try to be more clever */
               /* TODO: could we do something similar even if there are bilinear terms? */
               SCIP_Real f;

               f = SCIPfloor(scip, xval);
               rowcoef += consdata->quadsqrcoefs[j] * (2*f+1);
               bnd += consdata->quadsqrcoefs[j] * f * (f+1);
            }

            if( SCIPisInfinity(scip, ABS(rowcoef)) )
            {
               SCIPdebugMessage("skip linearization of square term in constraint %s because var %s is almost at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(x));
               SCIP_CALL( SCIPreleaseRow(scip, row) );
               return SCIP_OKAY;
            }
         }

         SCIP_CALL( SCIPaddVarToRow(scip, *row, x, rowcoef) );
      }

      for( j = 0; j < consdata->nbilinterms; ++j )
      { /* linearization of bilinear terms */
         coef = consdata->bilincoefs[j];
         x = consdata->bilinvars1[j];
         xval = SCIPgetSolVal(scip, sol, x);
         if( xval < SCIPvarGetLbLocal(x) )
            xval = SCIPvarGetLbLocal(x);
         else if( xval > SCIPvarGetUbLocal(x) )
            xval = SCIPvarGetUbLocal(x);
         if( SCIPisInfinity(scip, ABS(coef*xval)) )
         {
            SCIPdebugMessage("skip linearization of bilinear term in constraint %s because var %s is at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(x));
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            return SCIP_OKAY;
         }

         y = consdata->bilinvars2[j];
         yval = SCIPgetSolVal(scip, sol, y);
         if( yval < SCIPvarGetLbLocal(y) )
            yval = SCIPvarGetLbLocal(y);
         else if( yval > SCIPvarGetUbLocal(y) )
            yval = SCIPvarGetUbLocal(y);
         if( SCIPisInfinity(scip, ABS(coef*yval)) )
         {
            SCIPdebugMessage("skip linearization of bilinear term in constraint %s because var %s is at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(y));
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            return SCIP_OKAY;
         }

         SCIP_CALL( SCIPaddVarToRow(scip, *row, x, coef * yval) );
         SCIP_CALL( SCIPaddVarToRow(scip, *row, y, coef * xval) );
         bnd += coef * xval * yval;
      }
      
      if( SCIPisInfinity(scip, ABS(bnd)) )
      { /* seems to be a numerically bad cut */
         SCIPdebugMessage("skip cut for constraint %s because bound became very large: %g\n", SCIPconsGetName(cons), bnd);
         SCIP_CALL( SCIPreleaseRow(scip, row) );
         return SCIP_OKAY;
      }
     
      if( violbound == SCIP_BOUNDTYPE_LOWER )
      {
         SCIP_CALL( SCIPchgRowLhs(scip, *row, bnd) );
      }
      else
      {
         SCIP_CALL( SCIPchgRowRhs(scip, *row, bnd) );
      }
   }
   else
   { /* underestimate and linearize each term separately -> McCormick */
      for( j = 0; j < consdata->nquadvars; ++j )
      {
         rowcoef = consdata->quadlincoefs[j];

         x = consdata->quadvars[j];
         xval = SCIPgetSolVal(scip, sol, x);
         xlb = SCIPvarGetLbLocal(x);
         xub = SCIPvarGetUbLocal(x);
         if( xval < xlb )
            xval = xlb;
         else if( xval > xub )
            xval = xub;
         if( SCIPisInfinity(scip, ABS(xval)) )
         {
            SCIPdebugMessage("skip underestimator of square term in constraint %s because var %s is at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(x));
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            return SCIP_OKAY;
         }

         /* linearization of square term */
         coef = consdata->quadsqrcoefs[j];

         if( (violbound == SCIP_BOUNDTYPE_LOWER && coef <= 0) ||
             (violbound == SCIP_BOUNDTYPE_UPPER && coef >  0) )
         { /* convex -> linearize */
            if( SCIPvarGetType(x) == SCIP_VARTYPE_CONTINUOUS || SCIPisIntegral(scip, xval) )
            {
               rowcoef += 2*coef*xval;
               bnd += coef*xval*xval;
            }
            else
            { /* if variable is discrete but fractional, try to be more clever */
               SCIP_Real f = SCIPfloor(scip, xval);
               rowcoef += coef*(2*f+1);
               bnd += coef*f*(f+1);
            }
         }
         else
         { /* not convex -> secand approximation */
            if( SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub) )
            {
               SCIPdebugMessage("skip secand approx of square term in constraint %s because var %s is at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(x));
               SCIP_CALL( SCIPreleaseRow(scip, row) );
               return SCIP_OKAY;
            }

            rowcoef += coef * (xlb+xub);
            bnd += coef * xlb * xub;
         }

         if( SCIPisInfinity(scip, ABS(rowcoef)) )
         {
            SCIPdebugMessage("skip underestimator of square term in constraint %s because var %s is almost at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(x));
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            return SCIP_OKAY;
         }
         
         SCIP_CALL( SCIPaddVarToRow(scip, *row, x, rowcoef) );
      }

      for( j = 0; j < consdata->nbilinterms; ++j )
      {
         x = consdata->bilinvars1[j];
         xval = SCIPgetSolVal(scip, sol, x);
         xlb = SCIPvarGetLbLocal(x);
         xub = SCIPvarGetUbLocal(x);
         if( xval < xlb )
            xval = xlb;
         else if( xval > xub )
            xval = xub;
         if( SCIPisInfinity(scip, ABS(xval)) )
         {
            SCIPdebugMessage("skip underestimator of bilinear term in constraint %s because var %s is at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(x));
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            return SCIP_OKAY;
         }

         y = consdata->bilinvars2[j];
         yval = SCIPgetSolVal(scip, sol, y);
         ylb = SCIPvarGetLbLocal(y);
         yub = SCIPvarGetUbLocal(y);
         if( yval < ylb )
            yval = ylb;
         else if( yval > yub )
            yval = yub;
         if( SCIPisInfinity(scip, ABS(xval)) )
         {
            SCIPdebugMessage("skip underestimator of bilinear term in constraint %s because var %s is at infinity\n", SCIPconsGetName(cons), SCIPvarGetName(y));
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            return SCIP_OKAY;
         }

         coef = consdata->bilincoefs[j];
         if( violbound == SCIP_BOUNDTYPE_LOWER )
            coef = -coef;

         if( coef > 0.0 )
         {
            if( !SCIPisInfinity(scip, -xlb) && !SCIPisInfinity(scip, -ylb) &&
                (SCIPisInfinity(scip,  xub) ||  SCIPisInfinity(scip,  yub) ||
                  (xub-xlb)*yval + (yub-ylb)*xval <= xub*yub - xlb*ylb) )
            {
               xcoef = coef * ylb;
               ycoef = coef * xlb;
               bnd_  = coef * xlb * ylb;
            }
            else if( !SCIPisInfinity(scip, xub) && !SCIPisInfinity(scip, yub) )
            {
               xcoef = coef * yub;
               ycoef = coef * xub;
               bnd_  = coef * xub * yub;
            }
            else
            {
               SCIPdebugMessage("skip underestimator of bilinear term in constraint %s because var %s or %s is unbounded\n", SCIPconsGetName(cons), SCIPvarGetName(x), SCIPvarGetName(y));
               SCIP_CALL( SCIPreleaseRow(scip, row) );
               return SCIP_OKAY;
            }
         }
         else
         { /* coef < 0 */
            if( !SCIPisInfinity(scip,  xub) && !SCIPisInfinity(scip, -ylb) &&
                (SCIPisInfinity(scip, -xlb) ||  SCIPisInfinity(scip,  yub) ||
                  (xub-xlb)*yval - (yub-ylb)*xval <= xub*ylb - xlb*yub) )
            {
               xcoef = coef * ylb;
               ycoef = coef * xub;
               bnd_  = coef * xub * ylb;
            }
            else if( !SCIPisInfinity(scip, -xlb) && !SCIPisInfinity(scip, yub) )
            {
               xcoef = coef * yub;
               ycoef = coef * xlb;
               bnd_  = coef * xlb * yub;
            }
            else
            {
               SCIPdebugMessage("skip underestimator of bilinear term in constraint %s because var %s or %s is unbounded\n", SCIPconsGetName(cons), SCIPvarGetName(x), SCIPvarGetName(y));
               SCIP_CALL( SCIPreleaseRow(scip, row) );
               return SCIP_OKAY;
            }
         }
         if( violbound == SCIP_BOUNDTYPE_LOWER )
         {
            xcoef = -xcoef;
            ycoef = -ycoef;
            bnd_  = -bnd_;
         }
         
         if( SCIPisInfinity(scip, ABS(xcoef)) || SCIPisInfinity(scip, ABS(ycoef)) )
         {
            SCIPdebugMessage("skip underestimator of bilinear term in constraint %s because var %s or %s is almost unbounded\n", SCIPconsGetName(cons), SCIPvarGetName(x), SCIPvarGetName(y));
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            return SCIP_OKAY;
         }
         
         SCIP_CALL( SCIPaddVarToRow(scip, *row, x, xcoef) );
         SCIP_CALL( SCIPaddVarToRow(scip, *row, y, ycoef) );
         bnd += bnd_;
      }
      
      if( violbound == SCIP_BOUNDTYPE_LOWER )
      {
         SCIP_CALL( SCIPchgRowLhs(scip, *row, bnd) );
      }
      else
      {
         SCIP_CALL( SCIPchgRowRhs(scip, *row, bnd) );
      }
   }

   SCIPdebugMessage(" -> found cut rhs=%g, min=%f, max=%f range=%g\n",
       ABS(bnd),
       SCIPgetRowMinCoef(scip, *row), SCIPgetRowMaxCoef(scip, *row),
       SCIPgetRowMaxCoef(scip, *row)/SCIPgetRowMinCoef(scip, *row));

   if( SCIPisInfinity(scip, ABS(bnd)) )
   { /* seems to be a numerically bad cut */
      SCIPdebugMessage("skip cut for constraint %s because of very large left or right hand side: %g\n", SCIPconsGetName(cons), bnd);
      SCIP_CALL( SCIPreleaseRow(scip, row) );
      return SCIP_OKAY;
   }

   if( SCIPisGT(scip, SCIPgetRowMaxCoef(scip, *row)/SCIPgetRowMinCoef(scip, *row), maxrange) )
   { /* seems to be a numerically bad cut */
      SCIPdebugMessage("skip cut for constraint %s because of very large range: %g\n", SCIPconsGetName(cons), SCIPgetRowMaxCoef(scip, *row)/SCIPgetRowMinCoef(scip, *row));
      SCIP_CALL( SCIPreleaseRow(scip, row) );
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** tries to separate solution or LP solution by a linear cut
 * 
 *  assumes that constraint violations have been computed 
 */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< quadratic constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of constraints that seem to be useful */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   SCIP_RESULT*          result,             /**< result of separation */
   SCIP_Bool             addweakcuts         /**< should also weak (only slightly violated) cuts be added in a nonconvex constraint ? */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          efficacy;
   SCIP_BOUNDTYPE     violbound;
   int                c;
   SCIP_ROW*          row;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nusefulconss <= nconss);
   assert(result != NULL);
   
   *result = SCIP_FEASIBLE;
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( SCIPisFeasPositive(scip, consdata->lhsviol) || SCIPisFeasPositive(scip, consdata->rhsviol) )
      {
         /* we are not feasible anymore */
         if( *result == SCIP_FEASIBLE )
            *result = SCIP_DIDNOTFIND;

         violbound = SCIPisFeasPositive(scip, consdata->lhsviol) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER;

         /* generate cut */
         SCIP_CALL( generateCut(scip, conss[c], sol, violbound, &row, conshdlrdata->cutmaxrange) );
         
         /* if generation failed, then probably because of numerical issues;
          * if the constraint is convex and we are desperate to get a cut, then we can try again with a better chosen reference point */
         if( row == NULL && !addweakcuts &&
            ( (violbound == SCIP_BOUNDTYPE_UPPER && consdata->isconvex ) || /* convex  constraint, or */
              (violbound == SCIP_BOUNDTYPE_LOWER && consdata->isconcave) )  /* concave constraint */
           )
         {
            SCIP_CALL( generateCutCareful(scip, conss[c], sol, violbound, &row, conshdlrdata->cutmaxrange) );
         }
         
         if( row == NULL ) /* failed to generate cut */
            continue;

         efficacy = SCIPgetCutEfficacy(scip, sol, row);

         if( efficacy > conshdlrdata->mincutefficacy ||  /* ''strong'' cut */
             (SCIPisFeasPositive(scip, efficacy) &&       /* ''weak'' cut, use only if */ 
                (addweakcuts ||   /* flag is set, or */
                   (violbound == SCIP_BOUNDTYPE_UPPER && consdata->isconvex ) || /* convex  constraint, or */
                   (violbound == SCIP_BOUNDTYPE_LOWER && consdata->isconcave)    /* concave constraint */
                 )) )
         { /* cut cuts off solution */
            SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE /* forcecut */) );
            *result = SCIP_SEPARATED;
            SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
         }

         SCIP_CALL( SCIPreleaseRow (scip, &row) );
      }

      /* enforce only useful constraints
       * others are only checked and enforced if we are still feasible or have not found a separating cut yet
       */ 
      if( c >= nusefulconss && *result == SCIP_SEPARATED )
         break;
   }

   return SCIP_OKAY;
}

/** processes the event that a new primal solution has been found */
static
SCIP_DECL_EVENTEXEC(processNewSolutionEvent)
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS**    conss;
   int            nconss;
   SCIP_CONSDATA* consdata;
   int            c;
   SCIP_SOL*      sol;
   SCIP_ROW*      row = NULL;

   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   assert((SCIPeventGetType(event) | SCIP_EVENTTYPE_SOLFOUND) != 0);

   conshdlr = (SCIP_CONSHDLR*)eventdata;

   nconss = SCIPconshdlrGetNConss(conshdlr);

   if( nconss == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   sol = SCIPeventGetSol(event);
   assert(sol != NULL);

   /* we are only interested in solution coming from the NLP or RENSNL heuristic (is that good?) */
   if( SCIPsolGetHeur(sol) == NULL )
      return SCIP_OKAY;
   if( SCIPsolGetHeur(sol) != conshdlrdata->nlpheur && SCIPsolGetHeur(sol) != conshdlrdata->rensnlheur)
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(conshdlr);
   assert(conss != NULL);

   SCIPdebugMessage("catched new sol event %d from heur %p; have %d conss\n", SCIPeventGetType(event), (void*)SCIPsolGetHeur(sol), nconss);

   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsIsLocal(conss[c]) )
         continue;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->isconvex && !SCIPisInfinity(scip, consdata->rhs) )
      {
         SCIP_CALL( generateCut(scip, conss[c], sol, SCIP_BOUNDTYPE_UPPER, &row, conshdlrdata->cutmaxrange) );
      }
      else if( consdata->isconcave && !SCIPisInfinity(scip, -consdata->lhs) )
      {
         SCIP_CALL( generateCut(scip, conss[c], sol, SCIP_BOUNDTYPE_LOWER, &row, conshdlrdata->cutmaxrange) );
      }
      else
         continue;

      if( row == NULL )
         continue;

      assert(!SCIProwIsLocal(row));

      SCIP_CALL( SCIPaddPoolCut(scip, row) );
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   return SCIP_OKAY;
}

#ifdef USE_RELAXBRANCH
/** determines branching point for a variable */
static
SCIP_Real calculateBranchingPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             mindistbrpointtobound, /**< minimal (fractional) distance of branching point to variable bounds */
   SCIP_VAR*             var                 /**< variable to compute branching point for */
   )
{
   SCIP_Real branchpoint;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(var  != NULL);
   assert(mindistbrpointtobound >= 0.0);
   assert(mindistbrpointtobound <= 0.5);
     
   branchpoint = SCIPgetVarSol(scip, var);
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   if( SCIPisInfinity(scip, branchpoint) )
   {
      if( SCIPisPositive(scip, lb) )
         branchpoint = lb + 1000;
      else
         branchpoint = 0.0;
   }
   else if( SCIPisInfinity(scip, -branchpoint) )
   {
      if( SCIPisNegative(scip, ub) )
         branchpoint = ub - 1000;
      else
         branchpoint = 0.0;
   }

   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
   {
      if( !SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub) )
      {
         /* branch on value of LP solution
          * if it is too close to the bounds, move more into the middle of the interval */
         if( branchpoint < (1-mindistbrpointtobound) * lb + mindistbrpointtobound * ub )
            branchpoint = (1-mindistbrpointtobound) * lb + mindistbrpointtobound * ub;
         else if( branchpoint > mindistbrpointtobound * lb + (1-mindistbrpointtobound) * ub )
            branchpoint = mindistbrpointtobound * lb + (1-mindistbrpointtobound) * ub;

         /* for very tiny intervals we set it into the middle */
         if( SCIPisEQ(scip, lb/2, ub/2) )
            branchpoint = (lb+ub) * 0.5;
      }
      else if( !SCIPisLT(scip, lb, branchpoint) )
      {
         assert(SCIPisInfinity(scip, ub));
         branchpoint = lb + MAX(0.5*ABS(lb), 1000);
      }
      else if( !SCIPisGT(scip, ub, branchpoint) )
      {
         assert(SCIPisInfinity(scip, -lb));
         branchpoint = ub - MAX(0.5*ABS(ub), 1000);
      }
   }
   else
   {
      assert(ub - lb > 0.9);
      if( branchpoint > ub )
         branchpoint = ub - 0.5;
      else if( branchpoint < lb )
         branchpoint = lb + 0.5;
   }
   
   return branchpoint;
}
#endif

/** computes the infeasibilities of variables from the convexification gaps in the constraints and notifies the branching rule about them
 */
static
SCIP_RETCODE registerVariableInfeasibilities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check */
   int                   nconss,             /**< number of constraints to check */
   int*                  nnotify             /**< counter for number of notifications performed */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   int                c;
   int                j;
   SCIP_Real          xlb;
   SCIP_Real          xub;
   SCIP_Real          xval;
   SCIP_Real          ylb;
   SCIP_Real          yub;
   SCIP_Real          yval;
   SCIP_Real          gap;
   SCIP_Real          coef_;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
#ifdef WITH_CONSBRANCHNL
   assert(conshdlrdata->branchnl != NULL);
#endif
   
   *nnotify = 0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      SCIPdebugMessage("con %s violation: %g %g  convex: %u %u\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, consdata->isconvex, consdata->isconcave);
      
      if( !consdata->nquadvars )
         continue;
      
      if( (!SCIPisFeasPositive(scip, consdata->lhsviol) || consdata->isconcave) &&
          (!SCIPisFeasPositive(scip, consdata->rhsviol) || consdata->isconvex ) )
         continue;
      SCIPdebugMessage("con %s violation: %g %g  convex: %u %u\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, consdata->isconvex, consdata->isconcave);
      
      for( j = 0; j < consdata->nquadvars; ++j )
      { /* square terms */
         if( (SCIPisFeasPositive(scip, consdata->rhsviol) && consdata->quadsqrcoefs[j] < 0) ||
             (SCIPisFeasPositive(scip, consdata->lhsviol) && consdata->quadsqrcoefs[j] > 0) )
         {
            xlb = SCIPvarGetLbLocal(consdata->quadvars[j]);
            xub = SCIPvarGetUbLocal(consdata->quadvars[j]);
            if( SCIPisEQ(scip, xlb/2, xub/2) )
               continue;
            
            xval = SCIPgetSolVal(scip, NULL, consdata->quadvars[j]);
            
            if( SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub) )
               gap = SCIPinfinity(scip);
            else if( xval < xlb || xval > xub )
               continue;
            else
               gap = (xval-xlb)*(xub-xval)/(1+2*ABS(xval));
            assert(!SCIPisNegative(scip, gap));
#ifdef WITH_CONSBRANCHNL
            SCIP_CALL( SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(scip, conshdlrdata->branchnl, consdata->quadvars[j], MAX(gap, 0.), SCIPinfinity(scip)) );
#else
#ifdef USE_RELAXBRANCH
            SCIP_CALL( SCIPaddRelaxBranchCand(scip, consdata->quadvars[j], MAX(gap, 0.0), calculateBranchingPoint(scip, conshdlrdata->mindistbrpointtobound, consdata->quadvars[j])) );
#else
            SCIP_CALL( updateVarInfeasibility(scip, conshdlr, consdata->quadvars[j], MAX(gap, 0.0)) );
#endif
#endif
            ++*nnotify;
         }
      }

      for( j = 0; j < consdata->nbilinterms; ++j )
      { /* bilinear terms */
         xlb = SCIPvarGetLbLocal(consdata->bilinvars1[j]);
         xub = SCIPvarGetUbLocal(consdata->bilinvars1[j]);
         if( SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub) )
         {
#ifdef WITH_CONSBRANCHNL
            SCIP_CALL(  SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(scip, conshdlrdata->branchnl, consdata->bilinvars1[j], SCIPinfinity(scip), SCIPinfinity(scip)) );
#else
#ifdef USE_RELAXBRANCH
            SCIP_CALL( SCIPaddRelaxBranchCand(scip, consdata->bilinvars1[j], SCIPinfinity(scip), calculateBranchingPoint(scip, conshdlrdata->mindistbrpointtobound, consdata->bilinvars1[j])) );
#else
            SCIP_CALL( updateVarInfeasibility(scip, conshdlr, consdata->bilinvars1[j], SCIPinfinity(scip)) );
#endif
#endif
            ++*nnotify;
            continue;
         }

         ylb = SCIPvarGetLbLocal(consdata->bilinvars2[j]);
         yub = SCIPvarGetUbLocal(consdata->bilinvars2[j]);
         if( SCIPisInfinity(scip, -ylb) || SCIPisInfinity(scip, yub) )
         {
#ifdef WITH_CONSBRANCHNL
            SCIP_CALL(  SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(scip, conshdlrdata->branchnl, consdata->bilinvars2[j], SCIPinfinity(scip), SCIPinfinity(scip)) );
#else
#ifdef USE_RELAXBRANCH
            SCIP_CALL( SCIPaddRelaxBranchCand(scip, consdata->bilinvars2[j], SCIPinfinity(scip), calculateBranchingPoint(scip, conshdlrdata->mindistbrpointtobound, consdata->bilinvars2[j])) );
#else
            SCIP_CALL(  updateVarInfeasibility(scip, conshdlr, consdata->bilinvars2[j], SCIPinfinity(scip)) );
#endif
#endif
            ++*nnotify;
            continue;
         }

         xval = SCIPgetSolVal(scip, NULL, consdata->bilinvars1[j]);
         if( xval < xlb )
            xval = xlb;
         else if( xval > xub )
            xval = xub;
         
         yval = SCIPgetSolVal(scip, NULL, consdata->bilinvars2[j]);
         if( yval < ylb )
            yval = ylb;
         else if( yval > yub )
            yval = yub;
         
         coef_ = SCIPisFeasPositive(scip, consdata->lhsviol) ? -consdata->bilincoefs[j] : consdata->bilincoefs[j];
         if( coef_ > 0.0 )
         {
            if( (xub-xlb)*yval + (yub-ylb)*xval <= xub*yub - xlb*ylb )
               gap = (xval*yval - xlb*yval - ylb*xval + xlb*ylb) / (1+sqrt(xval*xval + yval*yval));
            else
               gap = (xval*yval - xval*yub - yval*xub + xub*yub) / (1+sqrt(xval*xval + yval*yval));
         }
         else
         { /* coef_ < 0 */
            if( (xub-xlb)*yval - (yub-ylb)*xval <= xub*ylb - xlb*yub )
               gap = -(xval*yval - xval*ylb - yval*xub + xub*ylb) / (1+sqrt(xval*xval + yval*yval));
            else
               gap = -(xval*yval - xval*yub - yval*xlb + xlb*yub) / (1+sqrt(xval*xval + yval*yval));
         }
         
         assert(!SCIPisNegative(scip, gap));
         if( gap < 0.0 ) 
            gap = 0.0;
         
         if( !SCIPisEQ(scip, xlb/2, xub/2) )
         {
#ifdef WITH_CONSBRANCHNL
            SCIP_CALL(  SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(scip, conshdlrdata->branchnl, consdata->bilinvars1[j], gap, SCIPinfinity(scip)) );
#else
#ifdef USE_RELAXBRANCH
            SCIP_CALL( SCIPaddRelaxBranchCand(scip, consdata->bilinvars1[j], gap, calculateBranchingPoint(scip, conshdlrdata->mindistbrpointtobound, consdata->bilinvars1[j])) );
#else
            SCIP_CALL( updateVarInfeasibility(scip, conshdlr, consdata->bilinvars1[j], gap) );
#endif
#endif
            ++*nnotify;
         }
         if( !SCIPisEQ(scip, ylb/2, yub/2) )
         {
#ifdef WITH_CONSBRANCHNL
            SCIP_CALL(  SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(scip, conshdlrdata->branchnl, consdata->bilinvars2[j], gap, SCIPinfinity(scip)) );
#else
#ifdef USE_RELAXBRANCH
            SCIP_CALL( SCIPaddRelaxBranchCand(scip, consdata->bilinvars2[j], gap, calculateBranchingPoint(scip, conshdlrdata->mindistbrpointtobound, consdata->bilinvars2[j])) );
#else
            SCIP_CALL( updateVarInfeasibility(scip, conshdlr, consdata->bilinvars2[j], gap) );
#endif
#endif
            ++*nnotify;
         }
      }
   }

   return SCIP_OKAY;
}

/** registers a variable from a violated constraint as branching candidate that has a large absolute value in the LP relaxation */
static
SCIP_RETCODE registerLargeLPValueVariableForBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_VAR**            brvar               /**< buffer to store branching variable */
   )
{
   SCIP_CONSDATA*      consdata;
   SCIP_Real           val;
   SCIP_Real           brvarval;
   int                 i;
   int                 c;
   
   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);
   
   *brvar = NULL;
   brvarval = 0.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      
      if( !SCIPisFeasPositive(scip, consdata->lhsviol) && !SCIPisFeasPositive(scip, consdata->rhsviol) )
         continue;
      
      for( i = 0; i < consdata->nquadvars; ++i )
      {
         /* do not propose fixed variables */
         if( SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->quadvars[i]), SCIPvarGetUbLocal(consdata->quadvars[i])) )
            continue;
         val = SCIPgetSolVal(scip, NULL, consdata->quadvars[i]);
         if( ABS(val) > brvarval )
         {
            brvarval = ABS(val);
            *brvar = consdata->quadvars[i];
         }
      }
   }
   
   if( *brvar != NULL )
   {
#ifdef WITH_CONSBRANCHNL
      SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata           != NULL);
      assert(conshdlrdata->branchnl != NULL);
      SCIP_CALL( SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(scip, conshdlrdata->branchnl, *brvar, brvarval, SCIPinfinity(scip)) );
#else
#ifdef USE_RELAXBRANCH
      SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata           != NULL);
      SCIP_CALL( SCIPaddRelaxBranchCand(scip, *brvar, brvarval, calculateBranchingPoint(scip, conshdlrdata->mindistbrpointtobound, *brvar)) );
#else
      SCIP_CALL( updateVarInfeasibility(scip, conshdlr, *brvar, brvarval) );
#endif
#endif
   }
   
   return SCIP_OKAY;
}

/** solves a linear equation \f$ b*x \in rhs \f$ and reduces bounds on x or deduces infeasibility if possible */
static
SCIP_RETCODE propagateBoundsLinearVar(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_CONS*            cons,               /**< constraint where we currently propagate */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_VAR*             var,                /**< variable which domain we might reduce */
   SCIP_Real             b,                  /**< linear coefficient of variable */
   SCIP_INTERVAL         rhs,                /**< right hand side */
   SCIP_RESULT*          result,             /**< result of propagation */
   int*                  nchgbds             /**< buffer where to add number of tightened bounds */
   )
{
   SCIP_Bool infeas;
   SCIP_Bool tightened;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(result != NULL);
   assert(nchgbds != NULL);

   SCIPintervalDivScalar(intervalinfty, &rhs, rhs, b);
   
   if( SCIPisInfinity(scip, SCIPintervalGetInf(rhs)) || SCIPisInfinity(scip, -SCIPintervalGetSup(rhs)) )
   { /* domain outside [-infty, +infty] -> declare node infeasible */
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }

   if( !SCIPisInfinity(scip, -SCIPintervalGetInf(rhs)) )
   {
      SCIP_CALL( SCIPtightenVarLb(scip, var, SCIPintervalGetInf(rhs), FALSE, &infeas, &tightened) );
      if( infeas )
      {
         SCIPdebugMessage("found %s infeasible due to domain propagation for linear variable %s\n", SCIPconsGetName(cons), SCIPvarGetName(var));
         *result = SCIP_CUTOFF;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         return SCIP_OKAY;
      }
      if( tightened )
      {
         SCIPdebugMessage("tightened lower bound of linear variable %s in constraint %s to %g\n", SCIPvarGetName(var), SCIPconsGetName(cons), SCIPvarGetLbLocal(var));
         ++*nchgbds;
         *result = SCIP_REDUCEDDOM;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }

   if( !SCIPisInfinity(scip, SCIPintervalGetSup(rhs)) )
   {
      SCIP_CALL( SCIPtightenVarUb(scip, var, SCIPintervalGetSup(rhs), FALSE, &infeas, &tightened) );
      if( infeas )
      {
         SCIPdebugMessage("found %s infeasible due to domain propagation for linear variable %s\n", SCIPconsGetName(cons), SCIPvarGetName(var));
         *result = SCIP_CUTOFF;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         return SCIP_OKAY;
      }
      if( tightened )
      {
         SCIPdebugMessage("tightened upper bound of linear variable %s in constraint %s to %g\n", SCIPvarGetName(var), SCIPconsGetName(cons), SCIPvarGetUbLocal(var));
         ++*nchgbds;
         *result = SCIP_REDUCEDDOM;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }

   return SCIP_OKAY;
}

/** solves a quadratic equation \f$ a x^2 + b x \in rhs \f$ (with b an interval) and reduces bounds on x or deduces infeasibility if possible.
 */
static
SCIP_RETCODE propagateBoundsQuadVar(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_CONS*            cons,               /**< constraint where we currently propagate */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_VAR*             var,                /**< variable which bounds with might tighten */
   SCIP_Real             a,                  /**< coefficient in square term */
   SCIP_INTERVAL         b,                  /**< coefficient in linear term */
   SCIP_INTERVAL         rhs,                /**< right hand side of quadratic equation */
   SCIP_RESULT*          result,             /**< result of propagation */
   int*                  nchgbds             /**< buffer where to add number of tightened bounds */
   )
{
   SCIP_INTERVAL newrange;
   SCIP_Bool     infeas;
   SCIP_Bool     tightened;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(result != NULL);
   assert(nchgbds != NULL);

   if( SCIPvarGetLbLocal(var) >= 0.0 )
   { /* need only positive solutions */
      SCIPintervalSolveUnivariateQuadExpressionPositive(intervalinfty, &newrange, a, b, rhs);
   }
   else if( SCIPvarGetUbLocal(var) <= 0.0 )
   { /* need only negative solutions */
      SCIP_INTERVAL tmp;
      SCIPintervalSetBounds(&tmp, -SCIPintervalGetSup(b), -SCIPintervalGetInf(b));
      SCIPintervalSolveUnivariateQuadExpressionPositive(intervalinfty, &tmp, a, tmp, rhs);
      if( SCIPintervalIsEmpty(tmp) )
      {
         SCIPdebugMessage("found %s infeasible due to domain propagation for quadratic variable %s\n", SCIPconsGetName(cons), SCIPvarGetName(var));
         *result = SCIP_CUTOFF;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         return SCIP_OKAY;
      }
      SCIPintervalSetBounds(&newrange, -SCIPintervalGetSup(tmp), -SCIPintervalGetInf(tmp));
   }
   else
   {
      SCIPintervalSolveUnivariateQuadExpression(intervalinfty, &newrange, a, b, rhs);
   }

   if( SCIPisInfinity(scip, SCIPintervalGetInf(newrange)) || SCIPisInfinity(scip, -SCIPintervalGetSup(newrange)) )
   { /* domain outside [-infty, +infty] -> declare node infeasible */
      SCIPdebugMessage("found %s infeasible because propagated domain of quadratic variable %s is outside of (-infty, +infty)\n", SCIPconsGetName(cons), SCIPvarGetName(var));
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }

   if( SCIPintervalIsEmpty(newrange) )
   {
      SCIPdebugMessage("found %s infeasible due to domain propagation for quadratic variable %s\n", SCIPconsGetName(cons), SCIPvarGetName(var));
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   if( !SCIPisInfinity(scip, -SCIPintervalGetInf(newrange)) )
   {
      SCIP_CALL( SCIPtightenVarLb(scip, var, SCIPintervalGetInf(newrange), FALSE, &infeas, &tightened) );
      if( infeas )
      {
         SCIPdebugMessage("found %s infeasible due to domain propagation for quadratic variable %s\n", SCIPconsGetName(cons), SCIPvarGetName(var));
         *result = SCIP_CUTOFF;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         return SCIP_OKAY;
      }
      if( tightened )
      {
         SCIPdebugMessage("tightened lower bound of quadratic variable %s in constraint %s to %g\n", SCIPvarGetName(var), SCIPconsGetName(cons), SCIPvarGetLbLocal(var));
         ++*nchgbds;
         *result = SCIP_REDUCEDDOM;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }

   if( !SCIPisInfinity(scip, SCIPintervalGetSup(newrange)) )
   {
      SCIP_CALL( SCIPtightenVarUb(scip, var, SCIPintervalGetSup(newrange), FALSE, &infeas, &tightened) );
      if( infeas )
      {
         SCIPdebugMessage("found %s infeasible due to domain propagation for quadratic variable %s\n", SCIPconsGetName(cons), SCIPvarGetName(var));
         *result = SCIP_CUTOFF;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         return SCIP_OKAY;
      }
      if( tightened )
      {
         SCIPdebugMessage("tightened upper bound of quadratic variable %s in constraint %s to %g -> %g\n", SCIPvarGetName(var), SCIPconsGetName(cons), SCIPintervalGetSup(newrange), SCIPvarGetUbLocal(var));
         ++*nchgbds;
         *result = SCIP_REDUCEDDOM;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }

   return SCIP_OKAY;
}

/** updates the ranges for linear variables in a constraint data;
 *  adds up all ranges;
 *  returns whether there is a variable which range is entire
 */
static
void propagateBoundsUpdateLinRange(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_INTERVAL*        linrangesum,        /**< for summing up ranges of linear terms */
   int*                  entirevaridx        /**< buffer to store index of single variable which domain is entire, 
                                              *   or -1 if there is none, or -2 if there are at least two */ 
   )
{  /*lint --e{666}*/
   SCIP_VAR* var;
   int       i;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(linrangesum != NULL);
   assert(entirevaridx != NULL);

   SCIPintervalSet(linrangesum,  0.0);
   *entirevaridx = -1;
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      if( SCIPintervalIsEmpty(consdata->linrange[i]) )
      {
         var = consdata->linvars[i];
         SCIPintervalSetBounds(&consdata->linrange[i],
            -infty2infty(SCIPinfinity(scip), intervalinfty, -SCIPvarGetLbLocal(var)),
             infty2infty(SCIPinfinity(scip), intervalinfty,  SCIPvarGetUbLocal(var)));
         SCIPintervalMulScalar(intervalinfty, &consdata->linrange[i], consdata->linrange[i], consdata->lincoefs[i]);
      }
#ifndef NDEBUG
      else
      {
         SCIP_INTERVAL tmp;
         var = consdata->linvars[i];
         SCIPintervalSetBounds(&tmp,
            -infty2infty(SCIPinfinity(scip), intervalinfty, -SCIPvarGetLbLocal(var)),
             infty2infty(SCIPinfinity(scip), intervalinfty,  SCIPvarGetUbLocal(var)));
         SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->lincoefs[i]);
         assert(SCIPintervalIsSubsetEQ(intervalinfty, tmp, consdata->linrange[i]));
      }
#endif

      if( SCIPintervalIsEntire(intervalinfty, consdata->linrange[i]) )
      {
         if( *entirevaridx >= 0 )
         {
            *entirevaridx = -2;
            return;
         }
         *entirevaridx = i;
      }
      SCIPintervalAdd(intervalinfty, linrangesum, *linrangesum, consdata->linrange[i]);
   }
}

/** updates the ranges for quadratic terms associated to each variable in a constraint data (a*x^2+b*x);
 *  adds up all ranges;
 *  returns whether there is a term which range is entire
 */
static
void propagateBoundsUpdateQuadRangeVar(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_INTERVAL*        quadrangesum,       /**< for summing up ranges of quadratic terms */
   int*                  entirevaridx        /**< buffer to store index of single variable which domain is entire, 
                                              *   or -1 if there is none, or -2 if there are at least two */ 
   )
{  /*lint --e{666}*/
   SCIP_VAR*     var;
   SCIP_INTERVAL tmp;
   int           i;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(quadrangesum != NULL);
   assert(entirevaridx != NULL);

   SCIPintervalSet(quadrangesum, 0.0);
   *entirevaridx = -1; /* index of quadratic variable which range is entire */
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( SCIPintervalIsEmpty(consdata->quadrangevar[i]) )
      {
         var = consdata->quadvars[i];
         SCIPintervalSetBounds(&consdata->quadrangevar[i], 
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var))), 
             infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var))));

         if( consdata->quadlincoefs[i] )
         {
            SCIPintervalSet(&tmp, consdata->quadlincoefs[i]);
            SCIPintervalQuad(intervalinfty, &consdata->quadrangevar[i], consdata->quadsqrcoefs[i], tmp, consdata->quadrangevar[i]);
         }
         else
         {
            SCIPintervalSquare(intervalinfty, &consdata->quadrangevar[i], consdata->quadrangevar[i]);
            assert(SCIPintervalGetInf(consdata->quadrangevar[i]) < intervalinfty);
            SCIPintervalMulScalar(intervalinfty, &consdata->quadrangevar[i], consdata->quadrangevar[i], consdata->quadsqrcoefs[i]);
         }
      }
#ifndef NDEBUG
      else
      {
         var = consdata->quadvars[i];
         SCIPintervalSetBounds(&tmp, 
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var))), 
             infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var))));

         if( consdata->quadlincoefs[i] )
         {
            SCIP_INTERVAL tmp2;
            SCIPintervalSet(&tmp2, consdata->quadlincoefs[i]);
            SCIPintervalQuad(intervalinfty, &tmp, consdata->quadsqrcoefs[i], tmp2, tmp);
         }
         else
         {
            SCIPintervalSquare(intervalinfty, &tmp, tmp);
            SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->quadsqrcoefs[i]);
         }
         assert(SCIPintervalIsSubsetEQ(intervalinfty, tmp, consdata->quadrangevar[i]));
      }
#endif

      if( SCIPintervalIsEntire(intervalinfty, consdata->quadrangevar[i]) )
      {
         if( *entirevaridx >= 0 )
         {
            *entirevaridx = -2;
            return; /* cannot reduce bounds on any variable if there are more than one variable with an entire range */
         }
         *entirevaridx = i;
      }
      
      SCIPintervalAdd(intervalinfty, quadrangesum, *quadrangesum, consdata->quadrangevar[i]);
   }
}

/** updates the ranges for bilinear terms in a constraint data;
 *  adds up all ranges
 */
static
void propagateBoundsUpdateBilinRange(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_INTERVAL*        bilinrangesum       /**< for summing up ranges of bilinear terms */
   )
{  /*lint --e{666}*/
   SCIP_VAR*      var;
   SCIP_INTERVAL  tmp;
   int            i;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(bilinrangesum != NULL);

   SCIPintervalSet(bilinrangesum, 0.);
   for( i = 0; i < consdata->nbilinterms; ++i )
   { /* check if a bilinrange need an update */
      if( SCIPintervalIsEmpty(consdata->bilinrange[i]) )
      {
         var = consdata->bilinvars1[i];
         SCIPintervalSetBounds(&consdata->bilinrange[i],
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var))), 
             infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var))));
         var = consdata->bilinvars2[i];
         SCIPintervalSetBounds(&tmp,
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var))),
             infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var))));
         SCIPintervalMul(intervalinfty, &consdata->bilinrange[i], consdata->bilinrange[i], tmp);
         assert(SCIPintervalGetInf(consdata->bilinrange[i]) <  intervalinfty);
         assert(SCIPintervalGetSup(consdata->bilinrange[i]) > -intervalinfty);
         SCIPintervalMulScalar(intervalinfty, &consdata->bilinrange[i], consdata->bilinrange[i], consdata->bilincoefs[i]);
      }
#ifndef NDEBUG
      else 
      {
         SCIP_INTERVAL tmp2;
         var = consdata->bilinvars1[i];
         SCIPintervalSetBounds(&tmp2, 
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var))),
             infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var))));
         var = consdata->bilinvars2[i];
         SCIPintervalSetBounds(&tmp,
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var))),
             infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var))));
         SCIPintervalMul(intervalinfty, &tmp, tmp2, tmp);
         SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilincoefs[i]);
         assert(SCIPintervalIsSubsetEQ(intervalinfty, tmp, consdata->bilinrange[i]));
      }
#endif
      SCIPintervalAdd(intervalinfty, bilinrangesum, *bilinrangesum, consdata->bilinrange[i]);
   }
}

/** propagates bounds on a quadratic constraint */
static
SCIP_RETCODE propagateBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_RESULT*          result,             /**< pointer to store the result of the propagation call */
   int*                  nchgbds             /**< buffer where to add the the number of changed bounds */
  )
{  /*lint --e{666}*/
   SCIP_CONSDATA*     consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_INTERVAL      consbounds;    /* lower and upper bounds of constraint */
   SCIP_Real          intervalinfty; /* infinity used for interval computation */  
   int                entirevaridx;  /* index of a variable which domain is entire */
   SCIP_INTERVAL      linrangesum;   /* range of linear part */
   SCIP_INTERVAL      quadrangesum;  /* sum of ranges of quadratic variable parts (ax^2+bx) */
   SCIP_INTERVAL      bilinrangesum; /* range of complete bilinear part */
   SCIP_VAR*          var;
   SCIP_Real          a;   /* quadratic coefficient of quadratic equation */
   SCIP_INTERVAL      b;   /* linear coefficient of quadratic equation */ 
   SCIP_INTERVAL      rhs; /* right hand side of quadratic equation */
   SCIP_INTERVAL      tmp;
   int                i;
   int                j;
   int                k;
   int                l;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(result != NULL);
   assert(nchgbds != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->linrange     || consdata->nlinvars    == 0);
   assert(consdata->quadrangevar || consdata->nquadvars   == 0);
   assert(consdata->bilinrange   || consdata->nbilinterms == 0);

   *result = SCIP_DIDNOTRUN;

   if( consdata->ispropagated )
      return SCIP_OKAY;

#ifndef NDEBUG
   /* assert that there are no variables that are fixed to -/+ infinity */
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      assert(!SCIPisInfinity(scip,  SCIPvarGetLbLocal(consdata->linvars[i])));
      assert(!SCIPisInfinity(scip, -SCIPvarGetUbLocal(consdata->linvars[i])));
   }

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      assert(!SCIPisInfinity(scip,  SCIPvarGetLbLocal(consdata->quadvars[i])));
      assert(!SCIPisInfinity(scip, -SCIPvarGetUbLocal(consdata->quadvars[i])));
   }
#endif

   consdata->ispropagated = TRUE;

   *result = SCIP_DIDNOTFIND;
   intervalinfty = 1000 * SCIPinfinity(scip) * SCIPinfinity(scip);

   SCIPintervalSetBounds(&consbounds,
      -infty2infty(SCIPinfinity(scip), intervalinfty, -consdata->lhs+SCIPfeastol(scip)/2.0),
       infty2infty(SCIPinfinity(scip), intervalinfty,  consdata->rhs+SCIPfeastol(scip)/2.0));
   
   propagateBoundsUpdateLinRange(scip, consdata, intervalinfty, &linrangesum, &entirevaridx);
   if( entirevaridx == -2 )
      return SCIP_OKAY; /* at least two variables that are completely unbounded; cannot propagate anything */

   if( SCIPintervalIsEmpty(consdata->quadrange) )
   { /* quadrange needs update */
      SCIPintervalSet(&consdata->quadrange, 0.0);
      SCIP_CALL( addQuadRange(scip, cons, intervalinfty, &consdata->quadrange, NULL) );
   }

   if( entirevaridx >= 0 )
   {
      assert(entirevaridx < consdata->nlinvars);
      if( SCIPintervalIsEntire(intervalinfty, consdata->quadrange) )
         return SCIP_OKAY; /* nothing we can do */

      var = consdata->linvars[entirevaridx];
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR);

      rhs = consbounds;
      for( i = 0; i < consdata->nlinvars; ++i )
         if( i != entirevaridx )
            SCIPintervalSub(intervalinfty, &rhs, rhs, consdata->linrange[i]);
      SCIPintervalSub(intervalinfty, &rhs, rhs, consdata->quadrange);

      SCIP_CALL( propagateBoundsLinearVar(scip, cons, intervalinfty, var, consdata->lincoefs[entirevaridx], rhs, result, nchgbds) );

      return SCIP_OKAY;
   }
   if( SCIPintervalIsEntire(intervalinfty, linrangesum) ) /* can still happen if two half-unbounded were added -> bad luck */
      return SCIP_OKAY;

   /* intersects linrangesum+quadrange with consbounds */
   SCIPintervalAdd(intervalinfty, &tmp, linrangesum, consdata->quadrange);
#if 1 /* hopefully ok now where bug in interval arith. fixed */
   SCIPintervalIntersect(&consbounds, consbounds, tmp);
   if( SCIPintervalIsEmpty(consbounds) )
   {
      SCIPdebugMessage("found %s infeasible due to forward propagation\n", SCIPconsGetName(cons));
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }
#else
   SCIPintervalIntersect(&tmp, consbounds, tmp);
   if( SCIPintervalIsEmpty(tmp) )
   { /* check again with slightly larger bounds: workaround for instance product */
      SCIP_INTERVAL tmp2;
      SCIPintervalSetBounds(&tmp2, -SCIPfeastol(scip)/2, SCIPfeastol(scip)/2);
      SCIPintervalAdd(intervalinfty, &tmp, linrangesum, consdata->quadrange);
      SCIPintervalAdd(intervalinfty, &tmp, tmp, tmp2);
      SCIPintervalIntersect(&consbounds, consbounds, tmp);
      if( SCIPintervalIsEmpty(consbounds) )
      {
         SCIPdebugMessage("found %s infeasible due to forward propagation\n", SCIPconsGetName(cons));
         *result = SCIP_CUTOFF;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         return SCIP_OKAY;
      }
      else
      {
         SCIPdebugMessage("node only feasible after widening bounds by %g\n", SCIPfeastol(scip));
      }
   }
   else
      consbounds = tmp;
#endif
   
   /* domain propagation for linear variables */
   for(i = 0; i < consdata->nlinvars; ++i)
   {
      var = consdata->linvars[i];
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR);
      assert(!SCIPisZero(scip, consdata->lincoefs[i]));

      SCIPintervalSub(intervalinfty, &rhs, consbounds, linrangesum);
      SCIPintervalUndoSub(intervalinfty, &rhs, rhs, consdata->linrange[i]);
      SCIPintervalSub(intervalinfty, &rhs, rhs, consdata->quadrange);

      SCIP_CALL( propagateBoundsLinearVar(scip, cons, intervalinfty, var, consdata->lincoefs[i], rhs, result, nchgbds) );
      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;
   }

   if( !consdata->nquadvars )
      return SCIP_OKAY;

   propagateBoundsUpdateQuadRangeVar(scip, consdata, intervalinfty, &quadrangesum, &entirevaridx);
   if( entirevaridx == -2 )
      return SCIP_OKAY; /* cannot reduce bounds on any variable if there are more than one variable which domain is entire */

   if( entirevaridx >= 0 )
   { /* there is exactly one quadratic variable which domain is entire; this one is the only chance where we can achieve a bound tightening then */ 
      var = consdata->quadvars[entirevaridx];
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR);

      SCIPintervalSub(intervalinfty, &rhs, consbounds, linrangesum);
      
      SCIPintervalSet(&tmp, 0.0);
      SCIP_CALL( addQuadRange(scip, cons, intervalinfty, &tmp, var) );
      SCIPintervalSub(intervalinfty, &rhs, rhs, tmp);

      /* add up coefficient interval in linear term of expression */
      SCIPintervalSet(&b, consdata->quadlincoefs[entirevaridx]);
      for( k = 0; k < consdata->nadjbilin[entirevaridx]; ++k )
      {
         l = consdata->adjbilin[entirevaridx][k];
         if( consdata->bilinvars1[l] == var )
         {
            assert(consdata->bilinvars2[l] != var);
            SCIPintervalSetBounds(&tmp, 
               -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvars2[l]), SCIPvarGetUbLocal(consdata->bilinvars2[l]))),
                infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvars2[l]), SCIPvarGetUbLocal(consdata->bilinvars2[l]))));
         }
         else
         {
            assert(consdata->bilinvars2[l] == var);
            SCIPintervalSetBounds(&tmp, 
               -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvars1[l]), SCIPvarGetUbLocal(consdata->bilinvars1[l]))),
                infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvars1[l]), SCIPvarGetUbLocal(consdata->bilinvars1[l]))));
         }
         SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilincoefs[l]);
         SCIPintervalAdd(intervalinfty, &b, b, tmp);
      }
      if( SCIPintervalIsEntire(SCIPinfinity(scip), b) )
         return SCIP_OKAY; /* no hope to reduce a bound */

      a = consdata->quadsqrcoefs[entirevaridx];

      SCIP_CALL( propagateBoundsQuadVar(scip, cons, intervalinfty, var, a, b, rhs, result, nchgbds) );    
      return SCIP_OKAY;
   }

   propagateBoundsUpdateBilinRange(scip, consdata, intervalinfty, &bilinrangesum);
   if( SCIPintervalIsEntire(intervalinfty, bilinrangesum) ) /* propagation on quad. vars makes no sense */
      return SCIP_OKAY;

   /* move everything into consbounds */
   SCIPintervalSub(intervalinfty, &consbounds, consbounds, linrangesum);
   if( conshdlrdata->fastpropagate && consdata->nquadvars > 2 )
   {
      SCIPintervalSub(intervalinfty, &consbounds, consbounds, bilinrangesum);
      SCIPintervalSub(intervalinfty, &consbounds, consbounds, quadrangesum);
   }

   if( SCIPintervalIsEntire(intervalinfty, consbounds) )
      return SCIP_OKAY;

   for( j = 0; j < consdata->nquadvars; ++j )
   {
      var = consdata->quadvars[j];
      assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR);

      /* skip fixed variables */
      if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         continue;

      /* setup equation a*x^2 + b*x \in rhs, x is var */
      a = consdata->quadsqrcoefs[j];

      SCIPintervalSet(&b, consdata->quadlincoefs[j]);
      if( conshdlrdata->fastpropagate )
      {
         if( consdata->nquadvars > 2 )
            SCIPintervalUndoSub(intervalinfty, &rhs, consbounds, consdata->quadrangevar[j]);
         else if( consdata->nquadvars == 2 )
         {
            if( SCIPintervalIsEmpty(consdata->quadrangevar[1-j]) )
            { /* this can happen if j==1 and we just improved the bound for j==0 */
               SCIPintervalSetBounds(&tmp,
                  -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->quadvars[1-j]), SCIPvarGetUbLocal(consdata->quadvars[1-j]))),
                   infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->quadvars[1-j]), SCIPvarGetUbLocal(consdata->quadvars[1-j]))));
               if( consdata->quadlincoefs[1-j] )
               {
                  SCIP_INTERVAL tmp2;
                  SCIPintervalSet(&tmp2, consdata->quadlincoefs[1-j]);
                  SCIPintervalQuad(intervalinfty, &tmp, consdata->quadsqrcoefs[1-j], tmp2, tmp);
               }
               else
               {
                  SCIPintervalSquare(intervalinfty, &tmp, tmp);
                  assert(SCIPintervalGetInf(tmp) < intervalinfty);
                  SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->quadsqrcoefs[1-j]);
               }
               SCIPintervalSub(intervalinfty, &rhs, consbounds, tmp);
            }
            else
            {
               SCIPintervalSub(intervalinfty, &rhs, consbounds, consdata->quadrangevar[1-j]);
            }
         }
         else
            rhs = consbounds;
         if( consdata->nadjbilin[j] > 0 )
         {
            /* we should not just put all bilinear terms into the right hand side, that would be fatal for equations like x*y \in [...]
             * since recomputing the best quad range is too expensive,
             * we undo all the substractions of bilinrange where this var is involved and setup an appropriate linear term (b)
             * 
             * this is also still expensive in case of many bilinear terms
             * thus, for now, we just pick some bilinear terms arbitrarily that are substracted from the rhs and leave the others on the rhs
             * TODO something better should be done, maybe sorting adjbilin with some criterion and take the first 6; or rewrite the whole propagation stuff
             */
#define MAXBILIN 6
            int step;
            if( consdata->nadjbilin[j] < MAXBILIN )
            {
               k = 0;
               step = 1;
            }
            else
            {
               k = SCIPgetDepth(scip) % MAXBILIN;
               step = consdata->nadjbilin[j] / MAXBILIN;
               assert(k < consdata->nadjbilin[j]);
               assert(step > 0);
            }
            for( ; k < consdata->nadjbilin[j]; k += step )
            {
               l = consdata->adjbilin[j][k];
               if( consdata->nquadvars > 2 )
               {
                  if( SCIPintervalIsEmpty(consdata->bilinrange[l]) )
                  { /* this might happen, if we just reduced a bound on a variable in this bilinear term
                     * however, we should not update bilinrange[l] here since it is not invalidated in future bound changes as long as the corresponding quadrangevar's are invalid
                     */
                     SCIP_INTERVAL tmp2;
                     SCIPintervalSetBounds(&tmp2,
                        -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvars1[l]), SCIPvarGetUbLocal(consdata->bilinvars1[l]))),
                         infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvars1[l]), SCIPvarGetUbLocal(consdata->bilinvars1[l]))));
                     SCIPintervalSetBounds(&tmp,
                        -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvars2[l]), SCIPvarGetUbLocal(consdata->bilinvars2[l]))),
                         infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvars2[l]), SCIPvarGetUbLocal(consdata->bilinvars2[l]))));
                     SCIPintervalMul(intervalinfty, &tmp, tmp, tmp2);
                     assert(SCIPintervalGetInf(tmp) <  intervalinfty);
                     assert(SCIPintervalGetSup(tmp) > -intervalinfty);
                     SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilincoefs[l]);
                     SCIPintervalUndoSub(intervalinfty, &rhs, rhs, tmp);
                  }
                  else
                     SCIPintervalUndoSub(intervalinfty, &rhs, rhs, consdata->bilinrange[l]);
               }

               if( consdata->bilinvars1[l] == var )
               {
                  assert(consdata->bilinvars2[l] != var);
                  SCIPintervalSetBounds(&tmp, 
                     -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvars2[l]), SCIPvarGetUbLocal(consdata->bilinvars2[l]))),
                      infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvars2[l]), SCIPvarGetUbLocal(consdata->bilinvars2[l]))));
               }
               else
               {
                  assert(consdata->bilinvars2[l] == var);
                  SCIPintervalSetBounds(&tmp, 
                     -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvars1[l]), SCIPvarGetUbLocal(consdata->bilinvars1[l]))),
                      infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvars1[l]), SCIPvarGetUbLocal(consdata->bilinvars1[l]))));
               }
               SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilincoefs[l]);
               SCIPintervalAdd(intervalinfty, &b, b, tmp);
            }
            if( SCIPintervalIsEntire(intervalinfty, b) )
               continue; /* no hope to reduce a bound */
         }
      }
      else
      {
         SCIPintervalSet(&tmp, 0.);
         SCIP_CALL( addQuadRange(scip, cons, intervalinfty, &tmp, var) ); /* add up everything in quad.part not belonging to var */
         SCIPintervalSub(intervalinfty, &rhs, consbounds, tmp);   /* put bounds - linrangesum - quadratic_except_var into rhs */
         for( k = 0; k < consdata->nadjbilin[j]; ++k )
         {
            l = consdata->adjbilin[j][k];
            if( consdata->bilinvars1[l] == var )
            {
               assert(consdata->bilinvars2[l] != var);
               SCIPintervalSetBounds(&tmp,
                  -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvars2[l]), SCIPvarGetUbLocal(consdata->bilinvars2[l]))),
                   infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvars2[l]), SCIPvarGetUbLocal(consdata->bilinvars2[l]))));
            }
            else
            {
               assert(consdata->bilinvars2[l] == var);
               SCIPintervalSetBounds(&tmp, 
                  -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinvars1[l]), SCIPvarGetUbLocal(consdata->bilinvars1[l]))),
                   infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinvars1[l]), SCIPvarGetUbLocal(consdata->bilinvars1[l]))));
            }
            SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilincoefs[l]);
            SCIPintervalAdd(intervalinfty, &b, b, tmp);
         }
         if( SCIPintervalIsEntire(intervalinfty, b) )
            continue; /* no hope to reduce a bound */
      }

      if( SCIPintervalIsEntire(intervalinfty, rhs) )
         continue;

      SCIP_CALL( propagateBoundsQuadVar(scip, cons, intervalinfty, var, a, b, rhs, result, nchgbds) );
      if( *result == SCIP_CUTOFF )
         break;
   }

   return SCIP_OKAY;
}

/** NLPI initialization method of constraint handler
 * 
 *  The constraint handler should create an NLPI representation of the constraints in the provided NLPI.
 */
SCIP_RETCODE SCIPconsInitNlpiQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for quadratic constraints */
   SCIP_NLPI*            nlpi,               /**< NLPI where to add constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_CONS**           conss,              /**< quadratic constraints */
   SCIP_HASHMAP*         scipvar2nlpvar      /**< mapping from SCIP variables to variable indices in NLPI */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real*     lhs;
   SCIP_Real*     rhs;
   int*           nlininds;
   int**          lininds;
   SCIP_Real**    linvals;
   int*           nquadrows;
   int**          quadrowidx;
   int**          quadoffset;
   int**          quadindex;
   SCIP_Real**    quadcoefs;
   int            quadnnz;
   SCIP_HASHMAP*  var2rowidx;
   SCIP_VAR*      othervar;
   int            i;
   int            j;
   int            k;
   int            l;
   int            lincnt;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(nlpi != NULL);

   if( nconss == 0 )
      return SCIP_OKAY;

   assert(conss != NULL);
   assert(scipvar2nlpvar != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nconss) );

   SCIP_CALL( SCIPallocBufferArray(scip, &nlininds, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds,  nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals,  nconss) );

   SCIP_CALL( SCIPallocBufferArray(scip, &nquadrows, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadrowidx, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadoffset, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadindex, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadcoefs, nconss) );

   for( i = 0; i < nconss; ++i )
   {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      
      /* skip local constraints; TODO do not add empty constraints to NLP */
      if( SCIPconsIsLocal(conss[i]) )
      {
         lhs[i] = -SCIPinfinity(scip);
         rhs[i] =  SCIPinfinity(scip);
         nlininds[i] = 0;
         nquadrows[i] = 0;
         continue;
      }

      lhs[i] = consdata->lhs;
      rhs[i] = consdata->rhs;
      
      /* count nonzeros in quadratic part */
      nlininds[i] = consdata->nlinvars;
      quadnnz = consdata->nbilinterms;
      for( j = 0; j < consdata->nquadvars; ++j )
      {
         if( consdata->quadsqrcoefs[j] )
            ++quadnnz;
         if( consdata->quadlincoefs[j] != 0.0 )
            ++nlininds[i];
      }

      if( nlininds[i] )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &lininds[i], nlininds[i]) );
         SCIP_CALL( SCIPallocBufferArray(scip, &linvals[i], nlininds[i]) );
      }
      else
      {
         lininds[i] = NULL;
         linvals[i] = NULL;
      }

      for( j = 0; j < consdata->nlinvars; ++j )
      {
         linvals[i][j] = consdata->lincoefs[j];
         assert(SCIPhashmapExists(scipvar2nlpvar, consdata->linvars[j]));
         lininds[i][j] = (int) (size_t) SCIPhashmapGetImage(scipvar2nlpvar, consdata->linvars[j]);
      }

      if( quadnnz == 0 )
      {
         nquadrows[i] = 0;
         quadrowidx[i] = NULL;
         quadoffset[i] = NULL;
         quadindex[i] = NULL;
         quadcoefs[i] = NULL;
         continue;
      }

      nquadrows[i] = consdata->nquadvars;
      SCIP_CALL( SCIPallocBufferArray(scip, &quadrowidx[i], consdata->nquadvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadoffset[i], consdata->nquadvars+1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadindex[i], quadnnz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadcoefs[i], quadnnz) );

      SCIP_CALL( SCIPhashmapCreate(&var2rowidx, SCIPblkmem(scip), consdata->nquadvars) );
      k = 0;
      lincnt = consdata->nlinvars;
      for( j = 0; j < consdata->nquadvars; ++j )
      {
         assert(SCIPhashmapExists(scipvar2nlpvar, consdata->quadvars[j]));
         if( consdata->quadlincoefs[j] != 0.0 )
         {
            lininds[i][lincnt] = (int)(size_t)SCIPhashmapGetImage(scipvar2nlpvar, consdata->quadvars[j]);
            linvals[i][lincnt] = consdata->quadlincoefs[j];
            ++lincnt;
         }

         assert( !SCIPhashmapExists(var2rowidx, consdata->quadvars[j]) );
         quadrowidx[i][j] = (int)(size_t)SCIPhashmapGetImage(scipvar2nlpvar, consdata->quadvars[j]);
         SCIP_CALL( SCIPhashmapInsert(var2rowidx, consdata->quadvars[j], (void*)(size_t)j) );

         quadoffset[i][j] = k;
         if( consdata->quadsqrcoefs[j] != 0.0 )
         {
            assert(k < quadnnz);
            quadindex[i][k] = j;
            quadcoefs[i][k] = consdata->quadsqrcoefs[j];
            ++k;
         }

         for( l = 0; l < consdata->nadjbilin[j]; ++l )
         {
            othervar = consdata->bilinvars1[consdata->adjbilin[j][l]];
            if( othervar == consdata->quadvars[j] )
               othervar = consdata->bilinvars2[consdata->adjbilin[j][l]];
            assert(othervar != consdata->quadvars[j]);

            if( SCIPhashmapExists(var2rowidx, othervar) ) /* processed the other var already, so now its time to add the corresponding bilin term */
            {
               assert(k < quadnnz);
               quadindex[i][k] = (int) (size_t) SCIPhashmapGetImage(var2rowidx, othervar);
               quadcoefs[i][k] = consdata->bilincoefs[consdata->adjbilin[j][l]];
               ++k;
            }
            /* otherwise we leave this bilinear term for later */
         }
      }
      quadoffset[i][consdata->nquadvars] = k;
      assert(k == quadnnz);
      SCIPhashmapFree(&var2rowidx);
   }

   SCIP_CALL( SCIPnlpiAddConstraints(scip, nlpi, nconss,
      lhs, rhs,
      nlininds, lininds, linvals,
      nquadrows, quadrowidx, quadoffset, quadindex, quadcoefs,
      NULL, NULL, NULL) );

   for( i = nconss-1; i >= 0; --i )
   {
      SCIPfreeBufferArrayNull(scip, &quadrowidx[i]);
      SCIPfreeBufferArrayNull(scip, &quadoffset[i]);
      SCIPfreeBufferArrayNull(scip, &quadindex[i]);
      SCIPfreeBufferArrayNull(scip, &quadcoefs[i]);
      SCIPfreeBufferArrayNull(scip, &lininds[i]);
      SCIPfreeBufferArrayNull(scip, &linvals[i]);
   }

   SCIPfreeBufferArray(scip, &quadcoefs);
   SCIPfreeBufferArray(scip, &quadindex);
   SCIPfreeBufferArray(scip, &quadoffset);
   SCIPfreeBufferArray(scip, &quadrowidx);
   SCIPfreeBufferArray(scip, &nquadrows);

   SCIPfreeBufferArray(scip, &nlininds);
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &linvals);

   SCIPfreeBufferArray(scip, &rhs);
   SCIPfreeBufferArray(scip, &lhs);

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int                i;
   
   assert(scip     != NULL);
   assert(conshdlr != NULL);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   for( i = 0; i < conshdlrdata->nquadconsupgrades; ++i )
   {
      assert(conshdlrdata->quadconsupgrades[i] != NULL);
      SCIPfreeMemory(scip, &conshdlrdata->quadconsupgrades[i]);
   }
   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->quadconsupgrades);
   
   SCIPfreeMemory(scip, &conshdlrdata);
   
   return SCIP_OKAY;
}

/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitQuadratic)
{  /*lint --e{715} */
   SCIP_CONSHDLRDATA* conshdlrdata;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

#ifdef WITH_CONSBRANCHNL
   conshdlrdata->branchnl = SCIPfindConshdlr(scip, "branchnonlinear");
   if( conshdlrdata->branchnl == NULL && nconss > 0 )
   {
      SCIPerrorMessage("cannot find constraint handler for branching on nonlinear variables");
      return SCIP_PLUGINNOTFOUND;
   }
#else
#ifndef USE_RELAXBRANCH
   /* TODO: what is a good estimate for the hashmap size? should the constraint handler notify about the number of potential candidates? */
   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->branchcand, SCIPblkmem(scip), SCIPgetNVars(scip) ? SCIPgetNVars(scip) : 1) );
#endif
#endif
   
#ifndef WITH_LAPACK
   if( nconss != 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Quadratic constraint handler does not have LAPACK for eigenvalue computation. Will assume that matrices (with size > 2x2) are indefinite.\n");
   }
#endif
#ifdef USECLOCK   
   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->clock1) );
   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->clock2) );
   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->clock3) );
#endif
   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitQuadratic)
{  /*lint --e{715} */
   SCIP_CONSHDLRDATA* conshdlrdata;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

#ifdef WITH_CONSBRANCHNL
   conshdlrdata->branchnl = NULL;
#else
#ifndef USE_RELAXBRANCH
   SCIP_CALL( clearBranchingCandidates(scip, conshdlr) );
   if( conshdlrdata->branchcand != NULL )
      SCIPhashmapFree(&conshdlrdata->branchcand);
#endif
#endif
   
#ifdef USECLOCK
   printf("clock1: %g\t clock2: %g\t clock3: %g\n", SCIPgetClockTime(scip, conshdlrdata->clock1), SCIPgetClockTime(scip, conshdlrdata->clock2), SCIPgetClockTime(scip, conshdlrdata->clock3));
   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->clock1) );
   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->clock2) );
   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->clock3) );
#endif
   return SCIP_OKAY;
}

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   *result = SCIP_FEASIBLE;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
   }

   return SCIP_OKAY;
}

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_HASHMAP*      terms;
   SCIP_Real          constant;
   SCIP_Bool          havechange;
   int                i;
   int                c;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;
   terms = NULL;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !consdata->isremovedfixings )
      {
         SCIP_CALL( presolveCreateQuadTerm(scip, &terms, &constant, &havechange, conss[c], conshdlrdata->replacesqrbinary) );

         if( havechange )
         {
            SCIPdebugMessage("exitpre found changed variable in constraint %s:\n", SCIPconsGetName(conss[c]));
            /* unlock all variables */
            for( i = 0; i < consdata->nlinvars; ++i )
            {
               SCIP_CALL( unlockLinearVariable(scip, conss[c], consdata->linvars[i], consdata->lincoefs[i]) );
            }
            for( i = 0; i < consdata->nquadvars; ++i )
            {
               SCIP_CALL( SCIPunlockVarCons(scip, consdata->quadvars[i], conss[c], TRUE, TRUE) );
            }
            SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );

            /* reset function */
            SCIP_CALL( consdataSetFunctionData(scip, consdata, terms, -1) );
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               consdata->lhs -= constant;
            if( !SCIPisInfinity(scip,  consdata->rhs) )
               consdata->rhs -= constant;

            /* lock all variables */
            for( i = 0; i < consdata->nlinvars; ++i )
            {
               SCIP_CALL( lockLinearVariable(scip, conss[c], consdata->linvars[i], consdata->lincoefs[i]) );
            }
            for( i = 0; i < consdata->nquadvars; ++i )
            {
               SCIP_CALL( SCIPlockVarCons(scip, consdata->quadvars[i], conss[c], TRUE, TRUE) );
            }
            SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
         }

         presolveQuadTermFree(scip, &terms);
      }

#ifndef NDEBUG
      for( i = 0; i < consdata->nlinvars; ++i )
         assert(SCIPvarIsActive(consdata->linvars[i]));

      for( i = 0; i < consdata->nquadvars; ++i )
         assert(SCIPvarIsActive(consdata->quadvars[i]));
#endif

      SCIP_CALL( boundUnboundedVars(scip, conss[c], conshdlrdata->defaultbound, NULL) );
   }

   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Bool          convex;
   int                c;
   
   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   convex = TRUE;
   
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      SCIP_CALL( checkCurvature(scip, conshdlr, conss[c]) );

      if( !SCIPisInfinity(scip,  consdata->rhs) && !consdata->isconvex )
      {
         SCIPdebugMessage("nonconvex because of upper bound in con %s\n", SCIPconsGetName(conss[c]));
         convex = FALSE;
      }
      if( !SCIPisInfinity(scip, -consdata->lhs) && !consdata->isconcave )
      {
         SCIPdebugMessage("nonconvex because of lower bound in con %s\n", SCIPconsGetName(conss[c]));
         convex = FALSE;
      }
   }
   
   if( nconss != 0 )
   {
       if( convex )
           SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "All quadratic constraints are convex\n");
       else
           SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "There are nonconvex quadratic constraints.\n");
   }

   conshdlrdata->nlpheur    = SCIPfindHeur(scip, "nlp");
   conshdlrdata->rensnlheur = SCIPfindHeur(scip, "rensnl");
   
   conshdlrdata->newsoleventfilterpos = -1;
   if( nconss != 0 && (conshdlrdata->nlpheur != NULL || conshdlrdata->rensnlheur != NULL) && conshdlrdata->linearizenlpsol )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_newsolution");
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, &conshdlrdata->newsoleventfilterpos) );
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int                c;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
   }

   if( conshdlrdata->newsoleventfilterpos >= 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;

      /* failing of the following events mean that new solution events should not have been catched */
      assert(conshdlrdata->nlpheur != NULL || conshdlrdata->rensnlheur != NULL);
      assert(conshdlrdata->linearizenlpsol);

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_newsolution");
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, conshdlrdata->newsoleventfilterpos) );
      conshdlrdata->newsoleventfilterpos = -1;
   }

   conshdlrdata->nlpheur    = NULL;
   conshdlrdata->rensnlheur = NULL;

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int                i;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );
   }

   for( i = 0; i < (*consdata)->nlinvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->linvars[i]) );
   }
   for( i = 0; i < (*consdata)->nquadvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->quadvars[i]) );
   }

   consdataFree(scip, consdata);

   assert(*consdata == NULL);

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransQuadratic)
{  
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
   int            i;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   SCIP_CALL( SCIPduplicateMemory(scip, &targetdata, sourcedata) );

   if( targetdata->nlinvars != 0 )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &targetdata->lincoefs, sourcedata->lincoefs, targetdata->nlinvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->linvars, targetdata->nlinvars) );
      SCIP_CALL( SCIPgetTransformedVars(scip, targetdata->nlinvars, sourcedata->linvars, targetdata->linvars) );
      for( i = 0; i < targetdata->nlinvars; ++i )
      {
         SCIP_CALL( SCIPcaptureVar(scip, targetdata->linvars[i]) );
      }
      SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->linrange, targetdata->nlinvars) );
   }
   else
   {
      targetdata->linvars = NULL;
      targetdata->lincoefs = NULL;
      targetdata->linrange = NULL;
   }
   targetdata->linbndchgeventdata  = NULL;
   targetdata->quadbndchgeventdata = NULL;

   if( targetdata->nquadvars != 0 )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &targetdata->quadlincoefs, sourcedata->quadlincoefs, targetdata->nquadvars) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &targetdata->quadsqrcoefs, sourcedata->quadsqrcoefs, targetdata->nquadvars) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &targetdata->nadjbilin,   sourcedata->nadjbilin,   targetdata->nquadvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->quadvars,  targetdata->nquadvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->adjbilin, targetdata->nquadvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->quadrangevar, targetdata->nquadvars) );
      SCIP_CALL( SCIPgetTransformedVars(scip, targetdata->nquadvars, sourcedata->quadvars, targetdata->quadvars) );
      for( i = 0; i < targetdata->nquadvars; ++i )
      {
         SCIP_CALL( SCIPcaptureVar(scip, targetdata->quadvars[i]) );
         if( targetdata->nadjbilin[i] != 0)
         {
            SCIP_CALL( SCIPduplicateMemoryArray(scip, &targetdata->adjbilin[i], sourcedata->adjbilin[i], sourcedata->nadjbilin[i]) );
         }
         else
            targetdata->adjbilin[i] = NULL;
      }

      if( targetdata->nbilinterms != 0)
      {
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &targetdata->bilincoefs, sourcedata->bilincoefs, targetdata->nbilinterms) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->bilinvars1, targetdata->nbilinterms) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->bilinvars2, targetdata->nbilinterms) );
         SCIP_CALL( SCIPgetTransformedVars(scip, targetdata->nbilinterms, sourcedata->bilinvars1, targetdata->bilinvars1) );
         SCIP_CALL( SCIPgetTransformedVars(scip, targetdata->nbilinterms, sourcedata->bilinvars2, targetdata->bilinvars2) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &targetdata->bilinrange, targetdata->nbilinterms) );
      }
      else
      {
         targetdata->bilincoefs = NULL;
         targetdata->bilinvars1 = NULL;
         targetdata->bilinvars2 = NULL;
         targetdata->bilinrange = NULL;
      }
   }
   else
   {
      targetdata->quadvars = NULL;
      targetdata->quadlincoefs = NULL;
      targetdata->quadsqrcoefs = NULL;
      targetdata->nadjbilin = NULL;
      targetdata->adjbilin = NULL;
      targetdata->bilincoefs = NULL;
      targetdata->bilinvars1 = NULL;
      targetdata->bilinvars2 = NULL;
      targetdata->bilinrange = NULL;
      targetdata->quadrangevar = NULL;
   }

   targetdata->isremovedfixings = FALSE;
   targetdata->ispropagated = FALSE;
   targetdata->ispresolved = FALSE;

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
      SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
      SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
      SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
      SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpQuadratic)
{  
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_ROW*          row;
   int                i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < nconss; ++i )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      if( !SCIPisInfinity(scip, -consdata->lhs) )
      {
         SCIP_CALL( generateCut(scip, conss[i], NULL, SCIP_BOUNDTYPE_LOWER, &row, conshdlrdata->cutmaxrange) );
         if( row != NULL )
         {
            SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE /* forcecut */) );
#ifdef SCIP_DEBUG
            SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
            SCIP_CALL( SCIPreleaseRow (scip, &row) );
         }
      }
      if( !SCIPisInfinity(scip, consdata->rhs) )
      {
         SCIP_CALL( generateCut(scip, conss[i], NULL, SCIP_BOUNDTYPE_UPPER, &row, conshdlrdata->cutmaxrange) );
         if( row != NULL )
         {
            SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE /* forcecut */) );
#ifdef SCIP_DEBUG
            SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
            SCIP_CALL( SCIPreleaseRow (scip, &row) );
         }
      }
   }

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpQuadratic)
{  
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_Bool          intervalfeas;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, conshdlrdata->doscaling, &maxviolcon) );
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   SCIP_CALL( areIntervalFeasible(scip, conss, nconss, maxviolcon, &intervalfeas) );
   if( !intervalfeas )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, result, TRUE) );
   if( *result == SCIP_SEPARATED )
      return SCIP_OKAY;

   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_Bool          intervalfeas;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(sol != NULL);
   assert(result != NULL);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( computeViolations(scip, conss, nconss, sol, conshdlrdata->doscaling, &maxviolcon) );
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   SCIP_CALL( areIntervalFeasible(scip, conss, nconss, maxviolcon, &intervalfeas) );
   if( !intervalfeas )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, result, TRUE) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_Bool          intervalfeas;
   SCIP_RESULT        separateresult;
   int                nnotify;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, conshdlrdata->doscaling, &maxviolcon) );
   if( maxviolcon == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }
   
   *result = SCIP_INFEASIBLE;

   SCIP_CALL( areIntervalFeasible(scip, conss, nconss, maxviolcon, &intervalfeas) );
   if( !intervalfeas )
   {
       *result = SCIP_CUTOFF;
       return SCIP_OKAY;
   }
   
   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, &separateresult, FALSE) );
   if( separateresult == SCIP_SEPARATED )
   {
      *result = SCIP_SEPARATED;
      return SCIP_OKAY;
   }

   /* we are not feasible, the whole node is not infeasible, and we cannot find a cut
    * -> collect variables for branching
    */
   
   SCIPdebugMessage("separation failed; max viol: %g+%g\n", SCIPconsGetData(maxviolcon)->lhsviol, SCIPconsGetData(maxviolcon)->rhsviol);

   SCIP_CALL( registerVariableInfeasibilities(scip, conshdlr, conss, nconss, &nnotify) );
   if( nnotify == 0 && !solinfeasible)
   { /* fallback: separation probably failed because of numerical difficulties with a convex constraint;
        if noone declared solution infeasible yet, try to resolve by branching */ 
      SCIP_VAR* brvar = NULL;
      SCIP_CALL( registerLargeLPValueVariableForBranching(scip, conshdlr, conss, nconss, &brvar) );
      if( brvar == NULL )
      {
         SCIPwarningMessage("Could not find any branching variable candidate. Cutting off node.\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      else
      {
         SCIPdebugMessage("Could not find any usual branching variable candidate. Proposed variable %s with LP value %g for branching.\n", SCIPvarGetName(brvar), SCIPgetSolVal(scip, NULL, brvar));
      }
   }

#ifndef WITH_CONSBRANCHNL
#ifndef USE_RELAXBRANCH
   SCIP_CALL( enforceByBranching(scip, conshdlr, result) );
   assert(*result == SCIP_BRANCHED);
#endif
#endif
   
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_CONSDATA*     consdata;
   SCIP_Bool          intervalfeas;
   int                c, i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, conshdlrdata->doscaling, &maxviolcon) );
   if( maxviolcon == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }
   
   *result = SCIP_INFEASIBLE;
      
   SCIP_CALL( areIntervalFeasible(scip, conss, nconss, maxviolcon, &intervalfeas) );
   if( !intervalfeas )
   {
       *result = SCIP_CUTOFF;
       return SCIP_OKAY;
   }
   
   /* we are not feasible and we cannot proof that the whole node is infeasible
    * -> collect variables for branching
    */

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      SCIPdebugMessage("con %s violation: %g %g  convex: %u %u\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, consdata->isconvex, consdata->isconcave);
      
      if( consdata->nquadvars == 0 )
         continue;
      
      if( !SCIPisFeasPositive(scip, consdata->lhsviol) && !SCIPisFeasPositive(scip, consdata->rhsviol) )
         continue;
      
      SCIPdebugMessage("con %s violation: %g %g\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol);
      
      for( i = 0; i < consdata->nquadvars; ++i )
         if( !SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->quadvars[i]), SCIPvarGetUbLocal(consdata->quadvars[i])) )
         {
#ifdef WITH_CONSBRANCHNL
            SCIP_CALL( SCIPconshdlrBranchNonlinearUpdateVarInfeasibility(scip, conshdlrdata->branchnl, consdata->quadvars[i], consdata->lhsviol + consdata->rhsviol, SCIPinfinity(scip)) );
#else
#ifdef USE_RELAXBRANCH
            SCIP_CALL( SCIPaddRelaxBranchCand(scip, consdata->quadvars[i], consdata->lhsviol + consdata->rhsviol, calculateBranchingPoint(scip, conshdlrdata->mindistbrpointtobound, consdata->quadvars[i])) );
#else
            SCIP_CALL( updateVarInfeasibility(scip, conshdlr, consdata->quadvars[i], consdata->lhsviol + consdata->rhsviol) );
#endif
#endif
         }
   }

#ifndef WITH_CONSBRANCHNL
#ifndef USE_RELAXBRANCH
   SCIP_CALL( enforceByBranching(scip, conshdlr, result) );
   assert(*result == SCIP_BRANCHED);
#endif
#endif
   
   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropQuadratic)
{
   SCIP_RESULT propresult;
   int         c;
   int         nchgbds = 0;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   for( c = 0; c < nconss && *result != SCIP_CUTOFF; ++c )
   {
      assert(conss != NULL);
      SCIP_CALL( propagateBounds(scip, conshdlr, conss[c], &propresult, &nchgbds) );
      if( propresult != SCIP_DIDNOTFIND && propresult != SCIP_DIDNOTRUN )
         *result = propresult;
      if( c >= nusefulconss && *result != SCIP_DIDNOTFIND )
         break;
   }

   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolQuadratic)
{  /*lint --e{715,788}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   int                c, i;
   SCIP_Bool          havechange;
   SCIP_Bool          keepnotpresolved;
   SCIP_HASHMAP*      terms;
   SCIP_Real          constant;
   SCIP_CONS*         upgdconslhs;
   SCIP_CONS*         upgdconsrhs;
   
   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);
   
   *result = SCIP_DIDNOTFIND;
   keepnotpresolved = FALSE;
   
   if( nrounds > 0 && nnewfixedvars == 0 && nnewupgdconss == 0 && nnewchgbds == 0 && nnewaggrvars == 0 && nnewchgvartypes == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   
   /** call upgrade method before aggregated variables are replaced or disaggregation takes place */
   havechange = FALSE;
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      for( i = 0; i < consdata->nlinvars; ++i )
      {
         SCIP_CALL( unlockLinearVariable(scip, conss[c], consdata->linvars[i], consdata->lincoefs[i]) );
      }
      SCIP_CALL( presolveUpgrade(scip, conshdlr, conss[c], &upgdconslhs, &upgdconsrhs) );
      if( upgdconslhs != NULL || upgdconsrhs != NULL )
      {
         /* add the upgraded constraint(s) to the problem */
         if( upgdconslhs != NULL && upgdconslhs != upgdconsrhs )
         {
            SCIP_CALL( SCIPaddCons(scip, upgdconslhs) );
            SCIP_CALL( SCIPreleaseCons(scip, &upgdconslhs) );
         }

         if( upgdconsrhs != NULL )
         {
            SCIP_CALL( SCIPaddCons(scip, upgdconsrhs) );
            SCIP_CALL( SCIPreleaseCons(scip, &upgdconsrhs) );
         }

         (*nupgdconss)++;
         *result = SCIP_SUCCESS;
         havechange = TRUE;

         /* delete upgraded constraint */
         if( SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, consdata->rhs) )
         {
            SCIPdebugMessage("delete constraint %s after upgrade\n", SCIPconsGetName(conss[c]));
            SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
            SCIP_CALL( SCIPdelCons(scip, conss[c]) );
            continue;
         }
      }
      for( i = 0; i < consdata->nlinvars; ++i )
      {
         SCIP_CALL( lockLinearVariable(scip, conss[c], consdata->linvars[i], consdata->lincoefs[i]) );
      }
   }
   if( havechange )
      return SCIP_OKAY;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      havechange = FALSE;
      if( !consdata->isremovedfixings || !consdata->ispresolved || (!consdata->ispropagated && conshdlrdata->replacebinaryprodlength) )
      {
         SCIPdebugMessage("process constraint %s\n", SCIPconsGetName(conss[c]));
         terms = NULL;
         SCIP_CALL( presolveCreateQuadTerm(scip, &terms, &constant, &havechange, conss[c], conshdlrdata->replacesqrbinary) );
#ifdef SCIP_DEBUG
         SCIPdebugMessage("%d; %g ", havechange, constant);
         presolveQuadTermPrint(scip, terms, NULL);
         SCIPdebugPrintf("\n");
#endif

         if( conshdlrdata->replacebinaryprodlength > 0 )
         {
            int nconsadded;

            nconsadded = 0;

            if( conshdlrdata->empathy4and == 2 )
            { /* user really likes AND, so give him */
               SCIP_CALL( presolveTryAddAND(scip, conss[c], terms, &nconsadded) );
               if( nconsadded != 0 )
               { /* does this count as an upgrade? */
                  *result = SCIP_SUCCESS;
                  havechange = TRUE;
               }
            }

            SCIP_CALL( presolveTryAddLinearReform(scip, conss[c], terms, &nconsadded, conshdlrdata->replacebinaryprodlength, conshdlrdata->empathy4and) );
            if( nconsadded != 0 )
            { /* does this count as an upgrade? */
               *result = SCIP_SUCCESS;
               havechange = TRUE;
            }
         }

         if( havechange || !consdata->ispresolved )
         {
            int ncomponents; /* number of quadratic blocks */

            if( SCIPhashmapIsEmpty(terms) )
            { /* all variables fixed or removed, constraint now constant */
               presolveQuadTermFree(scip, &terms);
               if( (!SCIPisInfinity(scip, -consdata->lhs) && SCIPisFeasLT(scip, consdata->lhs, constant    )) ||
                   (!SCIPisInfinity(scip,  consdata->rhs) && SCIPisFeasGT(scip, constant,     consdata->rhs)) )
               { /* constant is out of bounds */
                  SCIPdebugMessage("constraint %s is constant and infeasible\n", SCIPconsGetName(conss[c]));
                  SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
                  SCIP_CALL( SCIPdelCons(scip, conss[c]) );
                  *result = SCIP_CUTOFF;
                  break;
               }
               else
               { /* constant inside bounds, i.e., constraint always feasible */
                  SCIPdebugMessage("constraint %s is constant and feasible, removing\n", SCIPconsGetName(conss[c]));
                  SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) ); /* because only the presolved form is constant, in consdata we might still have variables */
                  SCIP_CALL( SCIPdelCons(scip, conss[c]) );
                  ++*ndelconss;
                  *result = SCIP_SUCCESS;
                  continue;
               }
            }
            
            SCIP_CALL( presolveFindComponents(scip, terms, &ncomponents) );
            assert(ncomponents >= 0);

            if( ncomponents == 0 )
            { /* upgrade to linear constraint */
               SCIP_HASHMAPLIST* list;
               SCIP_VAR*         var;
               PRESOLVEQUADTERM* term;
               SCIP_CONS*        lincon;
               
               SCIPdebugMessage("upgrade to linear constraint\n");
               SCIP_CALL( SCIPcreateConsLinear(scip, &lincon, SCIPconsGetName(conss[c]), 0, NULL, NULL,
                  SCIPisInfinity(scip, -consdata->lhs) ? consdata->lhs : consdata->lhs-constant,
                  SCIPisInfinity(scip,  consdata->rhs) ? consdata->rhs : consdata->rhs-constant,
                  SCIPconsIsInitial(conss[c]), SCIPconsIsSeparated(conss[c]), SCIPconsIsEnforced(conss[c]),
                  SCIPconsIsChecked(conss[c]), SCIPconsIsPropagated(conss[c]),  SCIPconsIsLocal(conss[c]),
                  SCIPconsIsModifiable(conss[c]), SCIPconsIsDynamic(conss[c]), SCIPconsIsRemovable(conss[c]),
                  SCIPconsIsStickingAtNode(conss[c])) );
               for( i = 0; i < SCIPhashmapGetNLists(terms); ++i )
                  for( list = SCIPhashmapGetList(terms, i); list != NULL; list = SCIPhashmapListGetNext(list) )
                  {
                     var = (SCIP_VAR*)SCIPhashmapListGetOrigin(list);
                     term = (PRESOLVEQUADTERM*)SCIPhashmapListGetImage(list);
                     assert(SCIPisZero(scip, term->sqrcoefs));
                     assert(term->bilin == NULL);
                     if( SCIPisZero(scip, term->lincoefs) )
                        continue;
                     SCIP_CALL( SCIPaddCoefLinear(scip, lincon, var, term->lincoefs) );
                  }

#ifdef SCIP_DEBUG
               SCIP_CALL( SCIPprintCons(scip, lincon, NULL) );
#endif
               SCIP_CALL( SCIPaddCons(scip, lincon) );
               SCIP_CALL( SCIPreleaseCons(scip, &lincon) );
               SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
               SCIP_CALL( SCIPdelCons(scip, conss[c]) );
               ++*nupgdconss;
               *result = SCIP_SUCCESS;
               presolveQuadTermFree(scip, &terms);
               continue;
            }

            if( havechange || (ncomponents > 1 && conshdlrdata->disaggregation) )
            {
               SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );

               /* do disaggregation only if nothing had changed,
                * so upgrade methods below have a chance to see the whole constraint
                */
               if( ncomponents == 1 || !conshdlrdata->disaggregation || havechange )
               {
                  /* unlock all variables */
                  for( i = 0; i < consdata->nlinvars; ++i )
                  {
                     SCIP_CALL( unlockLinearVariable(scip, conss[c], consdata->linvars[i], consdata->lincoefs[i]) );
                  }
                  for( i = 0; i < consdata->nquadvars; ++i )
                  {
                     SCIP_CALL( SCIPunlockVarCons(scip, consdata->quadvars[i], conss[c], TRUE, TRUE) );
                  }
                  
                  SCIP_CALL( consdataSetFunctionData(scip, consdata, terms, -1) );
                  if( !SCIPisInfinity(scip, -consdata->lhs) )
                     consdata->lhs -= constant;
                  if( !SCIPisInfinity(scip,  consdata->rhs) )
                     consdata->rhs -= constant;
                  
                  if( ncomponents > 1 && conshdlrdata->disaggregation )
                     keepnotpresolved = TRUE;
               }
               else
               {
#if 1
                  SCIP_CALL( presolveDisaggregate(scip, conshdlr, conss[c], terms, constant, ncomponents) );
                  SCIP_CALL( SCIPdelCons(scip, conss[c]) );
                  ++*nupgdconss;
                  *result = SCIP_SUCCESS;
                  presolveQuadTermFree(scip, &terms);
                  continue;
#else                  
                  /* unlock all variables */
                  for( i = 0; i < consdata->nlinvars; ++i )
                  {
                     SCIP_CALL( unlockLinearVariable(scip, conss[c], consdata->linvars[i], consdata->lincoefs[i]) );
                  }
                  for( i = 0; i < consdata->nquadvars; ++i )
                  {
                     SCIP_CALL( SCIPunlockVarCons(scip, consdata->quadvars[i], conss[c], TRUE, TRUE) );
                  }
                  SCIP_CALL( presolveDisaggregate(scip, conshdlr, conss[c], terms, constant, ncomponents) );
                  havechange = TRUE;
#endif
               }

               /* lock all variables */
               for( i = 0; i < consdata->nlinvars; ++i )
               {
                  SCIP_CALL( lockLinearVariable(scip, conss[c], consdata->linvars[i], consdata->lincoefs[i]) );
               }
               for( i = 0; i < consdata->nquadvars; ++i )
               {
                  SCIP_CALL( SCIPlockVarCons(scip, consdata->quadvars[i], conss[c], TRUE, TRUE) );
               }
               SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
            }
         }

         presolveQuadTermFree(scip, &terms);

         consdata->isremovedfixings = TRUE;
#if 0         
         /** let other constraint handlers try upgrading cons to a more specific quadratic constraint
          * since the locking of the linear variables changes if their is an upgrade, we unlock them all here (could be improved)
          */
         for( i = 0; i < consdata->nlinvars; ++i )
         {
            SCIP_CALL( unlockLinearVariable(scip, conss[c], consdata->linvars[i], consdata->lincoefs[i]) );
         }
         SCIP_CALL( presolveUpgrade(scip, conshdlr, conss[c], &upgdconslhs, &upgdconsrhs) );
         if( upgdconslhs != NULL || upgdconsrhs != NULL )
         {
            /* add the upgraded constraint(s) to the problem */
            if( upgdconslhs != NULL && upgdconslhs != upgdconsrhs )
            {
               SCIP_CALL( SCIPaddCons(scip, upgdconslhs) );
               SCIP_CALL( SCIPreleaseCons(scip, &upgdconslhs) );
            }

            if( upgdconsrhs != NULL )
            {
               SCIP_CALL( SCIPaddCons(scip, upgdconsrhs) );
               SCIP_CALL( SCIPreleaseCons(scip, &upgdconsrhs) );
            }

            (*nupgdconss)++;
            
            /* delete upgraded constraint */
            if( SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, consdata->rhs) )
            {
               SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[c]) );
               SCIP_CALL( SCIPdelCons(scip, conss[c]) );
               *result = SCIP_SUCCESS;
               continue;
            }
            else
               havechange = TRUE;
         }
         for( i = 0; i < consdata->nlinvars; ++i )
         {
            SCIP_CALL( lockLinearVariable(scip, conss[c], consdata->linvars[i], consdata->lincoefs[i]) );
         }
#endif
      }

      if( !consdata->ispropagated )
      { /* try domain propagation if there were bound changes or constraint has changed (in which case, catchVarEvents has set ispropagated to false) */
         SCIP_RESULT propresult;
         SCIP_CALL( propagateBounds(scip, conshdlr, conss[c], &propresult, nchgbds) );
         switch( propresult )
         {
            case SCIP_REDUCEDDOM:
               *result = SCIP_SUCCESS;
               break;
            case SCIP_CUTOFF:
               SCIPdebugMessage("propagation on constraint %s says problem is infeasible in presolve\n", SCIPconsGetName(conss[c]));
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            default:
               assert(propresult == SCIP_DIDNOTFIND || propresult == SCIP_DIDNOTRUN);
         }

         if( propresult != SCIP_REDUCEDDOM && !SCIPisInfinity(scip, conshdlrdata->defaultbound) )
         {
            if( nrounds > 0 )
            {
               int nboundchanges;

               nboundchanges = 0;

               SCIP_CALL( boundUnboundedVars(scip, conss[c], conshdlrdata->defaultbound, &nboundchanges) );
               if( nboundchanges != 0 )
               {
                  *nchgbds += nboundchanges;
                  *result   = SCIP_SUCCESS;
               }
            }
            else
            { /* wait for next round (or do in exitpre if no next round) */
               consdata->ispropagated = FALSE;
            }
         }
      }

      if( (nnewchgvartypes != 0 || havechange || !consdata->ispresolved)
          && (SCIPisEQ(scip, consdata->lhs, consdata->rhs) && SCIPisIntegral(scip, consdata->lhs)) )
      { /* check if we have a single linear continuous variable that we can make implicit integer 
        *  TODO allow for coefficient != +/-1 in front of linear var
        */
         int       ncontvar;
         SCIP_VAR* candidate;
         SCIP_Bool fail;

         fail = FALSE;
         candidate = NULL;
         ncontvar = 0;
         
         for( i = 0; !fail && i < consdata->nlinvars; ++i )
            if( !SCIPisIntegral(scip, consdata->lincoefs[i]) )
               fail = TRUE;
            else if( SCIPvarGetType(consdata->linvars[i]) == SCIP_VARTYPE_CONTINUOUS )
            {
               if( ncontvar > 0 ) /* now at 2nd continuous variable */
                  fail = TRUE;
               else if( SCIPisEQ(scip, ABS(consdata->lincoefs[i]), 1.0) )
                  candidate = consdata->linvars[i];
               ++ncontvar;
            }
         for( i = 0; !fail && i < consdata->nquadvars; ++i )
            fail = SCIPvarGetType(consdata->quadvars[i]) == SCIP_VARTYPE_CONTINUOUS ||
                  !SCIPisIntegral(scip, consdata->quadlincoefs[i]) ||
                  !SCIPisIntegral(scip, consdata->quadsqrcoefs[i]);
         for( i = 0; !fail && i < consdata->nbilinterms; ++i )
            fail = !SCIPisIntegral(scip, consdata->bilincoefs[i]);

         if( !fail && candidate != NULL )
         {
            SCIPdebugMessage("make variable %s implicit integer due to constraint %s\n", SCIPvarGetName(candidate), SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPchgVarType(scip, candidate, SCIP_VARTYPE_IMPLINT) );
            ++*nchgvartypes;
            *result = SCIP_SUCCESS;
         }
      }

      consdata->ispresolved = !keepnotpresolved;
   }

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropQuadratic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of quadratic constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropQuadratic NULL
#endif

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool      haslb;
   SCIP_Bool      hasub;
   int            i;
   
   assert(scip != NULL);
   assert(cons != NULL);
  
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   haslb = !SCIPisInfinity(scip, -consdata->lhs);
   hasub = !SCIPisInfinity(scip, consdata->rhs);

   for( i = 0; i < consdata->nlinvars; ++i )
   {
      if( consdata->lincoefs[i] > 0 )
      {
         if( haslb )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlockspos, nlocksneg) );
         }
         if( hasub )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlocksneg, nlockspos) );
         }
      }
      else
      {
         if( haslb )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlocksneg, nlockspos) );
         }
         if( hasub )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlockspos, nlocksneg) );
         }
      }
   }

   for( i = 0; i < consdata->nquadvars; ++i )
   { /* TODO try to be more clever */
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->quadvars[i], nlockspos+nlocksneg, nlockspos+nlocksneg) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveQuadratic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of quadratic constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveQuadratic NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveQuadratic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of quadratic constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveQuadratic NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableQuadratic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of quadratic constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableQuadratic NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableQuadratic)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of quadratic constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableQuadratic NULL
#endif


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int            j;
   
   assert(scip != NULL);
   assert(cons != NULL);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* print left hand side for ranged rows */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);

   /* print coefficients and variables */
   if( consdata->nlinvars == 0 && consdata->nquadvars == 0 )
      SCIPinfoMessage(scip, file, "0 ");
   else
   {
      for( j = 0; j < consdata->nlinvars; ++j )
      {
         SCIPinfoMessage(scip, file, "%+.15g<%s>[%c] ", consdata->lincoefs[j], SCIPvarGetName(consdata->linvars[j]),
            SCIPvarGetType(consdata->linvars[j]) == SCIP_VARTYPE_BINARY ? 'B' :
            SCIPvarGetType(consdata->linvars[j]) == SCIP_VARTYPE_INTEGER ? 'I' :
            SCIPvarGetType(consdata->linvars[j]) == SCIP_VARTYPE_IMPLINT ? 'I' : 'C');
      }

      for( j = 0; j < consdata->nquadvars; ++j )
      {
         if( consdata->quadlincoefs[j] != 0.0 )
         {
            SCIPinfoMessage(scip, file, "%+.15g<%s>[%c]", consdata->quadlincoefs[j], SCIPvarGetName(consdata->quadvars[j]),
               SCIPvarGetType(consdata->quadvars[j]) == SCIP_VARTYPE_BINARY ? 'B' :
               SCIPvarGetType(consdata->quadvars[j]) == SCIP_VARTYPE_INTEGER ? 'I' :
               SCIPvarGetType(consdata->quadvars[j]) == SCIP_VARTYPE_IMPLINT ? 'I' : 'C');
         }
         if( consdata->quadsqrcoefs[j] != 0.0 )
         {
            SCIPinfoMessage(scip, file, "%+.15g<%s>[%c]^2", consdata->quadsqrcoefs[j], SCIPvarGetName(consdata->quadvars[j]),
               SCIPvarGetType(consdata->quadvars[j]) == SCIP_VARTYPE_BINARY ? 'B' :
               SCIPvarGetType(consdata->quadvars[j]) == SCIP_VARTYPE_INTEGER ? 'I' :
               SCIPvarGetType(consdata->quadvars[j]) == SCIP_VARTYPE_IMPLINT ? 'I' : 'C',
               SCIPvarGetName(consdata->quadvars[j]),
               SCIPvarGetType(consdata->quadvars[j]) == SCIP_VARTYPE_BINARY ? 'B' :
               SCIPvarGetType(consdata->quadvars[j]) == SCIP_VARTYPE_INTEGER ? 'I' :
               SCIPvarGetType(consdata->quadvars[j]) == SCIP_VARTYPE_IMPLINT ? 'I' : 'C');
         }
      }

      for( j = 0; j < consdata->nbilinterms; ++j )
      {
         SCIPinfoMessage(scip, file, "%+.15g<%s>[%c]<%s>[%c]", consdata->bilincoefs[j],
            SCIPvarGetName(consdata->bilinvars1[j]),
            SCIPvarGetType(consdata->bilinvars1[j]) == SCIP_VARTYPE_BINARY ? 'B' :
            SCIPvarGetType(consdata->bilinvars1[j]) == SCIP_VARTYPE_INTEGER ? 'I' :
            SCIPvarGetType(consdata->bilinvars1[j]) == SCIP_VARTYPE_IMPLINT ? 'I' : 'C',
            SCIPvarGetName(consdata->bilinvars2[j]),
            SCIPvarGetType(consdata->bilinvars2[j]) == SCIP_VARTYPE_BINARY ? 'B' :
            SCIPvarGetType(consdata->bilinvars2[j]) == SCIP_VARTYPE_INTEGER ? 'I' :
            SCIPvarGetType(consdata->bilinvars2[j]) == SCIP_VARTYPE_IMPLINT ? 'I' : 'C');
      }
   }

   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
   {
      SCIPinfoMessage(scip, file, " == %.15g", consdata->rhs);
   }
   else if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      SCIPinfoMessage(scip, file, " <= %.15g", consdata->rhs);
   }
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      SCIPinfoMessage(scip, file, " >= %.15g", consdata->lhs);
   }
   else
   {
      SCIPinfoMessage(scip, file, " [free]");
   }

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          maxviol;
   int                c;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(sol != NULL);
   assert(result != NULL);
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;

   maxviol = 0.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      SCIP_CALL( computeViolation(scip, conss[c], sol, conshdlrdata->doscaling) );
      
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      
      if( SCIPisFeasPositive(scip, consdata->lhsviol) || SCIPisFeasPositive(scip, consdata->rhsviol) )
      {
         *result = SCIP_INFEASIBLE;
         if( printreason )
         {
            SCIPinfoMessage(scip, NULL, "quadratic constraint %s violated by %g+%g\n\t", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol);
            SCIP_CALL( consPrintQuadratic(scip, conshdlr, conss[c], NULL) );
         }
         if( conshdlrdata->nlpheur == NULL )
            return SCIP_OKAY;
         if( consdata->lhsviol > maxviol || consdata->rhsviol > maxviol )
            maxviol = consdata->lhsviol + consdata->rhsviol;
      }
   }
   
   if( *result == SCIP_INFEASIBLE && conshdlrdata->nlpheur )
   {
      SCIP_CALL( SCIPheurNlpUpdateStartpoint(scip, conshdlrdata->nlpheur, sol, maxviol) );
   }

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyQuadratic)
{  
   SCIP_CONSDATA* consdata;
   SCIP_VAR**     linvars;
   SCIP_VAR**     quadvars;
   SCIP_VAR**     bilinvars1;
   SCIP_VAR**     bilinvars2;
   int            i;
   int            j;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(sourcescip != NULL);
   assert(sourcecons != NULL);
   assert(varmap != NULL);
   assert(success != NULL);
   
   consdata = SCIPconsGetData(sourcecons);
   assert(consdata != NULL);
   
   *success = TRUE; /* think positive */

   linvars = NULL;
   quadvars = NULL;
   bilinvars1 = NULL;
   bilinvars2 = NULL;
   
   if( consdata->nlinvars != 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &linvars, consdata->nlinvars) );
      for( i = 0; i < consdata->nlinvars; ++i )
      {
         linvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, consdata->linvars[i]);
         if( linvars[i] == NULL )
         {
            *success = FALSE;
            break;
         }
      }
   }
   
   if( consdata->nbilinterms != 0 && *success )
   {
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &bilinvars1, consdata->nbilinterms) );
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &bilinvars2, consdata->nbilinterms) );
   }

   if( consdata->nquadvars != 0 && *success )
   {
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &quadvars, consdata->nquadvars) );
      for( i = 0; i < consdata->nquadvars; ++i )
      {
         quadvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, consdata->quadvars[i]);
         if( quadvars[i] == NULL )
         {
            *success = FALSE;
            break;
         }
         
         assert(consdata->nbilinterms != 0 || consdata->nadjbilin[i] == 0);
         
         for( j = 0; j < consdata->nadjbilin[i]; ++j )
         {
            assert(consdata->bilinvars1 != NULL);
            assert(consdata->bilinvars2 != NULL);
            assert(bilinvars1 != NULL);
            assert(bilinvars2 != NULL);
            if( consdata->bilinvars1[consdata->adjbilin[i][j]] == consdata->quadvars[i] )
            {
               assert(consdata->bilinvars2[consdata->adjbilin[i][j]] != consdata->quadvars[i]);
               bilinvars1[consdata->adjbilin[i][j]] = quadvars[i];
            }
            else
            {
               assert(consdata->bilinvars2[consdata->adjbilin[i][j]] == consdata->quadvars[i]);
               bilinvars2[consdata->adjbilin[i][j]] = quadvars[i];
            }
         }
      }
   }

   if(*success)
   {
      assert(stickingatnode == FALSE);
      SCIP_CALL( SCIPcreateConsQuadratic2(scip, cons, name ? name : SCIPconsGetName(sourcecons),
         consdata->nlinvars, linvars, consdata->lincoefs,
         consdata->nquadvars, quadvars, consdata->quadlincoefs, consdata->quadsqrcoefs,
         consdata->nadjbilin, consdata->adjbilin,
         consdata->nbilinterms, bilinvars1, bilinvars2, consdata->bilincoefs,
         consdata->lhs, consdata->rhs,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );
   }
   else
      *cons = NULL;
   

   SCIPfreeBufferArrayNull(sourcescip, &quadvars);
   SCIPfreeBufferArrayNull(sourcescip, &bilinvars2);
   SCIPfreeBufferArrayNull(sourcescip, &bilinvars1);
   SCIPfreeBufferArrayNull(sourcescip, &linvars);
   
   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
#if 1
static
SCIP_DECL_CONSPARSE(consParseQuadratic)
{  /*lint --e{715}*/
   SCIP_TOKENIZER tokenizer;
   SCIP_VAR** linvars;
   SCIP_Real* lincoefs;
   SCIP_VAR** bilinvars1;
   SCIP_VAR** bilinvars2;
   SCIP_Real* bilincoefs;
   SCIP_VAR*  var1;
   SCIP_VAR*  var2;

   SCIP_Real lhs;
   SCIP_Real rhs;

   SCIP_Bool havesign;
   SCIP_Bool havevalue;
   CIPSENSE sense;
   SCIP_Real coef;
   int coefsign;
   int lincoefssize;
   int nlincoefs;
   int bilincoefssize;
   int nbilincoefs;
   
   /* initialize tokenizer */
   tokenizer.strbuf = str;
   SCIP_CALL( SCIPallocBufferArray(scip, &tokenizer.token, SCIP_MAXSTRLEN) );
   tokenizer.token[0] = '\0';
   SCIP_CALL( SCIPallocBufferArray(scip, &tokenizer.tokenbuf, SCIP_MAXSTRLEN) );
   tokenizer.tokenbuf[0] = '\0';
   tokenizer.strpos = 0;
   
   /* initialize buffers for storing the coefficients */
   lincoefssize = 100;
   SCIP_CALL( SCIPallocBufferArray(scip, &linvars,  lincoefssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, lincoefssize) );
   bilincoefssize = 100;
   SCIP_CALL( SCIPallocBufferArray(scip, &bilinvars1, bilincoefssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bilinvars2, bilincoefssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bilincoefs, bilincoefssize) );
   
   /* set right hand and left side to their default values */
   lhs = -SCIPinfinity(scip);
   rhs =  SCIPinfinity(scip);

   (*success) = TRUE;

   /* read the coefficients */
   coefsign = +1;
   coef = 1.0;
   havesign  = FALSE;
   havevalue = FALSE;
   sense = CIP_SENSE_NOTHING;
   nlincoefs = 0;
   nbilincoefs = 0;
   var1 = NULL;
   var2 = NULL;

   SCIPdebugMessage("start parsing quadratic constraint expression\n");

   while( getNextToken(&tokenizer) && (*success) )
   {
      /* if the variable type is given ignore it */
      if( strncmp(tokenizer.token, "[B]", 3) == 0 && strlen(tokenizer.token) == 3 )
      {
         assert(var1 != NULL);
         assert(SCIPvarGetType(var1) == SCIP_VARTYPE_BINARY);
         SCIPdebugMessage("ignoring token <%s>\n", tokenizer.token);
         continue;
      }
      if( strncmp(tokenizer.token, "[I]", 3) == 0 && strlen(tokenizer.token) == 3 )
      {
         assert(var1 != NULL);
         assert(SCIPvarGetType(var1) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(var1) == SCIP_VARTYPE_IMPLINT);
         SCIPdebugMessage("ignoring token <%s>\n", tokenizer.token);
         continue;
      }
      if( strncmp(tokenizer.token, "[C]", 3) == 0 && strlen(tokenizer.token) == 3 )
      {
         assert(var1 != NULL);
         assert(SCIPvarGetType(var1) == SCIP_VARTYPE_CONTINUOUS);
         SCIPdebugMessage("ignoring token <%s>\n", tokenizer.token);
         continue;
      }
      
      if( *tokenizer.token == '>' && strlen(tokenizer.token) == 1 )
      {
         SCIPdebugMessage("ignoring token <%s>\n", tokenizer.token);
         continue;
      }
      if( *tokenizer.token == '<' && strlen(tokenizer.token) == 1 )
      {
         SCIPdebugMessage("ignoring token <%s>\n", tokenizer.token);
         continue;
      }
      
      /* check if we read a ^2 */
      if( strncmp(tokenizer.token, "^2", 2) == 0 && strlen(tokenizer.token) == 2 )
      {
         if( var1 == NULL )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error: expected variable before ^2 token\n");
            (*success) = FALSE;
            continue;
         }
         if( var2 != NULL )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error: unexpected ^2 token after <%s>*<%s>\n", SCIPvarGetName(var1), SCIPvarGetName(var2));
            (*success) = FALSE;
            continue;
         }
         SCIPdebugMessage("read ^2 symbol after variable <%s>\n", SCIPvarGetName(var1));
         var2 = var1;
         continue;
      }
      
      /* check if a term was finished */
      if( var1 != NULL )
      {
         CIPSENSE  nextsense = CIP_SENSE_NOTHING;
         int       nextcoefsign = +1;
         SCIP_Real nextcoef = 1.0;

         if( isSense(&tokenizer, &nextsense) || isSign(&tokenizer, &nextcoefsign) || isValue(scip, &tokenizer, &nextcoef) )
         {
            if( var2 == NULL )
            {
               /* insert the coefficient */
               SCIPdebugMessage("read coefficient: %+g<%s>\n", coefsign * coef, SCIPvarGetName(var1));
               if( !SCIPisZero(scip, coef) )
               {
                  /* resize the vars and coefs array if needed */
                  if( nlincoefs >= lincoefssize )
                  {
                     lincoefssize *= 2;
                     lincoefssize = MAX(lincoefssize, nlincoefs+1);
                     SCIP_CALL( SCIPreallocBufferArray(scip, &linvars,  lincoefssize) );
                     SCIP_CALL( SCIPreallocBufferArray(scip, &lincoefs, lincoefssize) );
                  }
                  assert(nlincoefs < lincoefssize);

                  /* add coefficient */
                  linvars [nlincoefs] = var1;
                  lincoefs[nlincoefs] = coefsign * coef;
                  nlincoefs++;
               }
            }
            else
            {
               /* insert the coefficient */
               SCIPdebugMessage("read coefficient: %+g<%s><%s>\n", coefsign * coef, SCIPvarGetName(var1), SCIPvarGetName(var2));
               if( !SCIPisZero(scip, coef) )
               {
                  /* resize the vars and coefs array if needed */
                  if( nbilincoefs >= bilincoefssize )
                  {
                     bilincoefssize *= 2;
                     bilincoefssize = MAX(bilincoefssize, nbilincoefs+1);
                     SCIP_CALL( SCIPreallocBufferArray(scip, &bilinvars1, bilincoefssize) );
                     SCIP_CALL( SCIPreallocBufferArray(scip, &bilinvars2, bilincoefssize) );
                     SCIP_CALL( SCIPreallocBufferArray(scip, &bilincoefs, bilincoefssize) );
                  }
                  assert(nbilincoefs < bilincoefssize);

                  /* add coefficient */
                  bilinvars1[nbilincoefs] = var1;
                  bilinvars2[nbilincoefs] = var2;
                  bilincoefs[nbilincoefs] = coefsign * coef;
                  nbilincoefs++;
               }
            }
            
            /* reset the flags and coefficient value for the next term */
            coefsign = +1;
            coef = 1.0;
            havesign = FALSE;
            havevalue = FALSE;
            var1 = NULL;
            var2 = NULL;
         }
      }
      
      /* check if we reached an equation sense */
      if( isSense(&tokenizer, &sense) )
      {
         if( nlincoefs == 0 && nbilincoefs == 0 && havevalue )
         {
            /* the constraint has no variables */
            switch(sense)
            {
            case CIP_SENSE_LE:
               lhs = coef;
               break;
            case CIP_SENSE_GE:
               rhs = coef;
               break;
            case CIP_SENSE_EQ:
               lhs = coef;
               rhs = coef;
               break;
            default:
               SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error: expected unknown sense <%d>)\n", sense);
               (*success) = FALSE;
            }
            sense = CIP_SENSE_NOTHING;
            havevalue = FALSE;
         }

         continue;
      }
      
      /* check if we read a sign */
      if( isSign(&tokenizer, &coefsign) )
      {
         SCIPdebugMessage("read coefficient sign: %+d\n", coefsign);
         havesign = TRUE;
         continue;
      }

      /* all but the first coefficient need a sign */
      if( (nlincoefs > 0 || nbilincoefs > 0) && !havesign && sense == CIP_SENSE_NOTHING )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error: expected sign ('+' or '-') or sense ('<' or '>')\n");
         (*success) = FALSE;
         break;
      }

      /* check if we read a value */
      if( isValue(scip, &tokenizer, &coef) )
      {
         SCIPdebugMessage("read coefficient value: <%g> with sign %+d\n", coef, coefsign);
         if( havevalue )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error: two consecutive values");
            (*success) = FALSE;
            break;
         }
         havevalue = TRUE;
         
         if( sense == CIP_SENSE_EQ ) 
         {
            lhs = coef;
            rhs = coef;
         }
         else if( sense == CIP_SENSE_LE ) 
            rhs = coef;
         if( sense == CIP_SENSE_GE ) 
            lhs = coef;
         
         continue;
      }
      
      /* the token is a variable name */
      
      /* if there is no first variable yet, get the corresponding variable */
      if( var1 == NULL )
      {
         var1 = SCIPfindVar(scip, tokenizer.token);

         if( var1 == NULL )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unknown variable <%s>", tokenizer.token);
            (*success) = FALSE;
         }
         
         continue;
      }
      
      /* if there is already a first and a second variable, then this is an error */
      if( var2 != NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "third variable variable <%s> in product", tokenizer.token);
         (*success) = FALSE;
         continue;
      }
      
      /* get the corresponding variable */
      var2 = SCIPfindVar(scip, tokenizer.token);

      if( var2 == NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unknown variable <%s>", tokenizer.token);
         (*success) = FALSE;
      }
   }
   
   if( var1 != NULL )
   { /* TODO may this happen for free constraints? */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unexpected end of string after term with variable <%s>", SCIPvarGetName(var1));
      (*success) = FALSE;
   }

   if( *success )
   {
      assert(stickingatnode == FALSE);
      SCIP_CALL( SCIPcreateConsQuadratic(scip, cons, name, nlincoefs, linvars, lincoefs, 
         nbilincoefs, bilinvars1, bilinvars2, bilincoefs, lhs, rhs, 
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );
   }

   SCIPfreeBufferArray(scip, &lincoefs);
   SCIPfreeBufferArray(scip, &linvars);
   SCIPfreeBufferArray(scip, &bilincoefs);
   SCIPfreeBufferArray(scip, &bilinvars1);
   SCIPfreeBufferArray(scip, &bilinvars2);
   tokenizer.token[0] = '\0';
   SCIPfreeBufferArray(scip, &tokenizer.tokenbuf);
   SCIPfreeBufferArray(scip, &tokenizer.token);
   
   return SCIP_OKAY;
}
#else
#define consParseQuadratic NULL
#endif


/*
 * constraint specific interface methods
 */

/** creates the handler for quadratic constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrQuadratic(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create quadratic constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   BMSclearMemory(conshdlrdata);

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         consFreeQuadratic, consInitQuadratic, consExitQuadratic,
         consInitpreQuadratic, consExitpreQuadratic, consInitsolQuadratic, consExitsolQuadratic,
         consDeleteQuadratic, consTransQuadratic, consInitlpQuadratic,
         consSepalpQuadratic, consSepasolQuadratic, consEnfolpQuadratic, consEnfopsQuadratic, consCheckQuadratic,
         consPropQuadratic, consPresolQuadratic, consRespropQuadratic, consLockQuadratic,
         consActiveQuadratic, consDeactiveQuadratic,
         consEnableQuadratic, consDisableQuadratic,
         consPrintQuadratic, consCopyQuadratic, consParseQuadratic,
         conshdlrdata) );

   /* add quadratic constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/replacesqrbinary",
         "whether a square of a binary variables should be replaced by the binary variable",
         &conshdlrdata->replacesqrbinary, FALSE, TRUE, NULL, NULL) );
/* 
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/replace_binaryprod",
         "whether products with a binary variable should be replaced by a new variable and an AND constraint
         (if possible) or a set of equivalent linear constraints", &conshdlrdata->replace_binaryprod, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/replace_binaryprod_forceAND",
         "whether products of binary variables should be replaced always by an AND constraint",
         &conshdlrdata->replace_binaryprod_forceAND, FALSE, FALSE, NULL, NULL) );
*/

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/replacebinaryprod",
         "max. length of linear term which when multiplied with a binary variables is replaced by an auxiliary variable and a linear reformulation (0 to turn off)",
         &conshdlrdata->replacebinaryprodlength, FALSE, INT_MAX, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/empathy4and",
         "empathy level for using the AND constraint handler: 0 always avoid using AND; 1 use AND sometimes; 2 use AND as often as possible",
         &conshdlrdata->empathy4and, FALSE, 0, 0, 2, NULL, NULL) );
   
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/disaggregate",
         "whether quadratic constraints consisting of several quadratic blocks should be disaggregated in several constraints",
         &conshdlrdata->disaggregation, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/minefficacy",
         "minimal efficacy for a cut to be added to the LP; overwrites separating/efficacy",
         &conshdlrdata->mincutefficacy, FALSE, 0.0001, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/scaling", 
         "whether a quadratic constraint should be scaled w.r.t. the current gradient norm when checking for feasibility",
         &conshdlrdata->doscaling, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/fastpropagate",
         "whether a propagation should be used that is faster in case of bilinear term, but also less efficient",
         &conshdlrdata->fastpropagate, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/defaultbound",
         "a default bound to impose on unbounded variables in quadratic terms (-defaultbound is used for missing lower bounds)",
         &conshdlrdata->defaultbound, TRUE, SCIPinfinity(scip), 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/cutmaxrange",
         "maximal range of a cut (maximal coefficient divided by minimal coefficient) in order to be added to LP relaxation",
         &conshdlrdata->cutmaxrange, FALSE, 1e+10, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/linearizenlpsol",
         "whether convex quadratic constraints should be linearized in a solution found by the NLP or RENSNL heuristic",
         &conshdlrdata->linearizenlpsol, FALSE, TRUE, NULL, NULL) );

#ifndef WITH_CONSBRANCHNL
#ifndef USE_RELAXBRANCH
   SCIP_CALL( SCIPaddCharParam(scip, "constraints/"CONSHDLR_NAME"/strategy",
         "strategy to use for selecting branching variable: b: rb-int-br, r: rb-int-br-rev, i: rb-inf",
         &conshdlrdata->strategy, FALSE, 'r', "bri", NULL, NULL) );
#endif

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/mindistbrpointtobound",
         "minimal fractional distance of branching point to variable bounds; a value of 0.5 leads to branching always in the middle of a bounded domain",
         &conshdlrdata->mindistbrpointtobound, FALSE, 0.2, 0.0001, 0.5, NULL, NULL) );
#endif
   
   SCIP_CALL( SCIPincludeEventhdlr(scip, CONSHDLR_NAME"_boundchange", "signals a bound change to a quadratic constraint",
      NULL, NULL, NULL, NULL, NULL, NULL, processVarEvent, NULL) );
   conshdlrdata->eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_boundchange");

   SCIP_CALL( SCIPincludeEventhdlr(scip, CONSHDLR_NAME"_newsolution", "handles the event that a new primal solution has been found",
      NULL, NULL, NULL, NULL, NULL, NULL, processNewSolutionEvent, NULL) );

   return SCIP_OKAY;
}

/** includes a quadratic constraint update method into the quadratic constraint handler */
SCIP_RETCODE SCIPincludeQuadconsUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_QUADCONSUPGD((*quadconsupgd)),  /**< method to call for upgrading quadratic constraint */
   int                   priority,           /**< priority of upgrading method */
   const char*           conshdlrname        /**< name of the constraint handler */
   )
{
   SCIP_CONSHDLR*        conshdlr;
   SCIP_CONSHDLRDATA*    conshdlrdata;
   SCIP_QUADCONSUPGRADE* quadconsupgrade;
   char                  paramname[SCIP_MAXSTRLEN];
   char                  paramdesc[SCIP_MAXSTRLEN];
   int                   i;
   
   assert(quadconsupgd != NULL);
   assert(conshdlrname != NULL );

   /* find the linear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("quadratic constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create a quadratic constraint upgrade data object */
   SCIP_CALL( SCIPallocMemory(scip, &quadconsupgrade) );
   quadconsupgrade->quadconsupgd = quadconsupgd;
   quadconsupgrade->priority     = priority;
   quadconsupgrade->active       = TRUE;

   /* insert quadratic constraint update method into constraint handler data */
   assert(conshdlrdata->nquadconsupgrades <= conshdlrdata->quadconsupgradessize);
   if( conshdlrdata->nquadconsupgrades+1 > conshdlrdata->quadconsupgradessize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, conshdlrdata->nquadconsupgrades+1);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->quadconsupgrades, newsize) );
      conshdlrdata->quadconsupgradessize = newsize;
   }
   assert(conshdlrdata->nquadconsupgrades+1 <= conshdlrdata->quadconsupgradessize);
   
   for( i = conshdlrdata->nquadconsupgrades; i > 0 && conshdlrdata->quadconsupgrades[i-1]->priority < quadconsupgrade->priority; --i )
      conshdlrdata->quadconsupgrades[i] = conshdlrdata->quadconsupgrades[i-1];
   assert(0 <= i && i <= conshdlrdata->nquadconsupgrades);
   conshdlrdata->quadconsupgrades[i] = quadconsupgrade;
   conshdlrdata->nquadconsupgrades++;

   /* adds parameter to turn on and off the upgrading step */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/"CONSHDLR_NAME"/upgrade/%s", conshdlrname);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "enable quadratic upgrading for constraint handler <%s>", conshdlrname);
   SCIP_CALL( SCIPaddBoolParam(scip,
         paramname, paramdesc,
         &quadconsupgrade->active, FALSE, TRUE, NULL, NULL) );
   
   return SCIP_OKAY;
}

/** Creates and captures a quadratic constraint.
 * 
 *  The constraint should be given in the form
 *  \f[
 *  \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m a_j y_jz_j \leq u,
 *  \f]
 *  where \f$x_i = y_j = z_k\f$ is possible.
 */
SCIP_RETCODE SCIPcreateConsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) */
   int                   nquadterms,         /**< number of quadratic terms (m) */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms (y_j) */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms (z_j) */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms (a_j) */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation (ell) */
   SCIP_Real             rhs,                /**< right hand side of quadratic equation (u) */
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
                                              *   are seperated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_HASHMAP*  terms;
   int i;
   
   assert( modifiable == FALSE ); /* we do not support column generation */

   /* find the quadratic constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("quadratic constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocMemory( scip, &consdata) );
   BMSclearMemory(consdata);

   consdata->lhs = lhs;
   consdata->rhs = rhs;

   SCIP_CALL( SCIPhashmapCreate(&terms, SCIPblkmem(scip), nlinvars + nquadterms) );
   for( i = 0; i < nlinvars; ++i )
   {
      if( SCIPisZero(scip, lincoefs[i]) )
         continue;
      SCIP_CALL( presolveQuadTermAdd(scip, terms, linvars[i], lincoefs[i], 0., NULL, 0.) );
   }
   for( i = 0; i < nquadterms; ++i )
   {
      if( SCIPisZero(scip, quadcoefs[i]) )
         continue;
      assert(SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED || SCIPvarGetStatus(quadvars1[i]) != SCIP_VARSTATUS_MULTAGGR);
      assert(SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED || SCIPvarGetStatus(quadvars2[i]) != SCIP_VARSTATUS_MULTAGGR);
      if( quadvars1[i] == quadvars2[i] )
      {
         SCIP_CALL( presolveQuadTermAdd(scip, terms, quadvars1[i], 0., quadcoefs[i], NULL, 0.) );
      }
      else
      {
         SCIP_CALL( presolveQuadTermAdd(scip, terms, quadvars1[i], 0., 0., quadvars2[i], quadcoefs[i]) );
         SCIP_CALL( presolveQuadTermAdd(scip, terms, quadvars2[i], 0., 0., quadvars1[i], quadcoefs[i]) );
      }
   }
   SCIP_CALL( consdataSetFunctionData(scip, consdata, terms, -1) );
   presolveQuadTermFree(scip, &terms);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, FALSE) );

   if( SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVING )
   {
      SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);
      
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, *cons) );

      if( SCIPgetStage(scip) > SCIP_STAGE_INITSOLVE )
      {
         SCIP_CALL( checkCurvature(scip, conshdlr, *cons) );
         SCIPdebugMessage("new quadratic constraint %s is %sconvex and %sconcave\n", name, consdata->isconvex ? "" : "not ", consdata->isconcave ? "" : "not ");
      }
   }

   return SCIP_OKAY;
}

/** Creates and captures a quadratic constraint.
 * 
 * The constraint should be given in the form
 * \f[
 * \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m (a_j y_j^2 + b_j y_j) + \sum_{k=1}^p c_kv_kw_k \leq u.
 * \f]
 */
SCIP_RETCODE SCIPcreateConsQuadratic2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) */ 
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) */ 
   int                   nquadvars,          /**< number of quadratic terms (m) */
   SCIP_VAR**            quadvars,           /**< array with variables in quadratic terms (y_j) */
   SCIP_Real*            quadlincoefs,       /**< array with linear coefficients of quadratic variables (b_j) */
   SCIP_Real*            quadsqrcoefs,       /**< array with coefficients of square terms of quadratic variables (a_j) */
   int*                  nadjbilin,          /**< array with numbers of bilinear terms where the variables are involved */
   int**                 adjbilin,           /**< array with indices of bilinear terms in which variable is involved */
   int                   nbilinterms,        /**< number of bilinear terms (p) */
   SCIP_VAR**            bilinvars1,         /**< array with first variables in bilinear terms (v_k) */
   SCIP_VAR**            bilinvars2,         /**< array with second variables in bilinear terms (w_k) */
   SCIP_Real*            bilincoefs,         /**< array with coefficients of bilinear terms (c_k) */
   SCIP_Real             lhs,                /**< constraint left hand side (ell) */
   SCIP_Real             rhs,                /**< constraint right hand side (u) */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint dynamic? */
   SCIP_Bool             removable           /**< should the constraint be removed from the LP due to aging or cleanup? */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   int i;
   
   assert( modifiable == FALSE ); /* we do not support column generation */
   assert(nlinvars == 0 || (linvars != NULL && lincoefs != NULL));
   assert(nquadvars == 0 || (quadvars != NULL && quadlincoefs != NULL && quadsqrcoefs != NULL));
   assert(nbilinterms == 0 || (bilinvars1 != NULL && bilinvars2 != NULL && bilincoefs != NULL && nquadvars > 0));
   
   /* find the quadratic constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("quadratic constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocMemory( scip, &consdata) );
   BMSclearMemory(consdata);

   consdata->lhs = lhs;
   consdata->rhs = rhs;

   if( nlinvars != 0 )
   {
      assert(linvars != NULL);
      consdata->nlinvars = nlinvars;
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->linvars,   linvars,   nlinvars) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->lincoefs, lincoefs, nlinvars) );
      for( i = 0; i < nlinvars; ++i )
      {
         SCIP_CALL( SCIPcaptureVar(scip, linvars[i]) );
      }
      SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->linrange, nlinvars) );
   }

   if( nquadvars != 0 )
   {
      assert(quadvars != NULL);
      consdata->nquadvars = nquadvars;
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->quadvars, quadvars, nquadvars) );
      if( quadlincoefs != NULL )
      {
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->quadlincoefs, quadlincoefs, nquadvars) );
      }
      else
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->quadlincoefs, nquadvars) );
         BMSclearMemoryArray(consdata->quadlincoefs, nquadvars);
      }
      if( quadsqrcoefs != NULL )
      {
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->quadsqrcoefs, quadsqrcoefs, nquadvars) );
      }
      else
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->quadsqrcoefs, nquadvars) );
         BMSclearMemoryArray(consdata->quadsqrcoefs, nquadvars);
      }
      if( nadjbilin != NULL )
      {
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->nadjbilin,   nadjbilin,   nquadvars) );
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->adjbilin,     adjbilin,   nquadvars) );
      }
      else
      {
         assert(nbilinterms == 0);
         SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->nadjbilin, nquadvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->adjbilin,  nquadvars) );
         BMSclearMemoryArray(consdata->nadjbilin, nquadvars);
         BMSclearMemoryArray(consdata->adjbilin,  nquadvars);
      }
      for( i = 0; i < nquadvars; ++i )
      {
         SCIP_CALL( SCIPcaptureVar(scip, quadvars[i]) );
         if( consdata->nadjbilin[i] != 0 )
         {
            assert(adjbilin[i] != NULL);
            assert(nadjbilin   != NULL);
            SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->adjbilin[i], adjbilin[i], nadjbilin[i]) );
         }
         else
         {
            assert(consdata->adjbilin[i] == NULL);
         }
      }
      SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->quadrangevar, nquadvars) );
   }

   if( nbilinterms != 0 )
   {
      consdata->nbilinterms = nbilinterms;
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->bilincoefs, bilincoefs, nbilinterms) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->bilinvars1,  bilinvars1,  nbilinterms) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &consdata->bilinvars2,  bilinvars2,  nbilinterms) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->bilinrange, nbilinterms) );
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
      local, modifiable, dynamic, removable, FALSE) );

   if( SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVING )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlrdata= SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);
      
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, *cons) );

      if( SCIPgetStage(scip) > SCIP_STAGE_INITSOLVE )
      {
         SCIP_CALL( checkCurvature(scip, conshdlr, *cons) );
         SCIPdebugMessage("new quadratic constraint %s is %sconvex and %sconcave\n", name, consdata->isconvex ? "" : "not ", consdata->isconcave ? "" : "not ");
      }
   }

   return SCIP_OKAY;
}

/** Gets the number of variables in the linear term of a quadratic constraint.
 */
int SCIPgetNLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->nlinvars;
}

/** Gets the variables in the linear part of a quadratic constraint.
 *  Length is given by SCIPgetNLinearVarsQuadratic.
 */
SCIP_VAR** SCIPgetLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->linvars;
}

/** Gets the number of variables in the quadratic term of a quadratic constraint.
 */
int SCIPgetNQuadVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->nquadvars;
}

/** Gets the variables in the quadratic part of a quadratic constraint.
 *  Length is given by SCIPgetNQuadVarsQuadratic.
 */
SCIP_VAR** SCIPgetQuadVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->quadvars;
}

/** Gets the coefficients in the linear part of a quadratic constraint.
 *  Length is given by SCIPgetNLinearVarsQuadratic.
 */
SCIP_Real* SCIPgetCoefsLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->lincoefs;
}

/** Gets the linear coefficients in the quadratic part of a quadratic constraint.
 *  Length is given by SCIPgetNQuadVarsQuadratic.
 */
SCIP_Real* SCIPgetLinearCoefsQuadVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->quadlincoefs;
}

/** Gets the square coefficients in the quadratic part of a quadratic constraint.
 *  Length is given by SCIPgetNQuadVarsQuadratic.
 */
SCIP_Real* SCIPgetSqrCoefsQuadVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->quadsqrcoefs;
}

/** Gets the number of bilinear terms in a quadratic constraint.
 */
int SCIPgetNBilinTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->nbilinterms;
}

/** Gets the first variables in the bilinear terms in a quadratic constraint.
 *  Length is given by SCIPgetNBilinTermQuadratic.
 */
SCIP_VAR** SCIPgetBilinVars1Quadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->bilinvars1;
}

/** Gets the second variables in the bilinear terms in a quadratic constraint.
 *  Length is given by SCIPgetNBilinTermQuadratic.
 */
SCIP_VAR** SCIPgetBilinVars2Quadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->bilinvars2;
}

/** Gets the coefficients of the bilinear terms in a quadratic constraint.
 *  Length is given by SCIPgetNBilinTermQuadratic.
 */
SCIP_Real* SCIPgetBilinCoefsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->bilincoefs;
}

/** Gets for each quadratic variable the number of bilinear terms in which the variable is involved in a quadratic constraint.
 *  Length is given by SCIPgetNQuadVarsQuadratic
 */
int* SCIPgetNAdjBilinQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->nadjbilin;
}

/** Gets for each quadratic variable the indices of bilinear terms in which the variable is involved in a quadratic constraint.
 *  Length is given by SCIPgetNQuadVarsQuadratic, length of each entry is given by SCIPgetNAdjBilinQuadratic.
 */
int** SCIPgetAdjBilinQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->adjbilin;
}

/** Gets the left hand side of a quadratic constraint.
 */
SCIP_Real SCIPgetLhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->lhs;
}

/** Gets the right hand side of a quadratic constraint.
 */
SCIP_Real SCIPgetRhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->rhs;
}

/** Indicates whether the quadratic function of a quadratic constraint is (known to be) convex.
 */
SCIP_Bool SCIPisConvexQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->isconvex;
}

/** Indicates whether the quadratic function of a quadratic constraint is (known to be) concave.
 */
SCIP_Bool SCIPisConcaveQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   
   return SCIPconsGetData(cons)->isconcave;
}

/** Computes the violation of a constraint by a solution */
SCIP_RETCODE SCIPgetViolationQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution which violation to calculate, or NULL for LP solution */
   SCIP_Real*            violation           /**< buffer to store violation of constraint */
   )
{
   SCIP_CONSDATA* consdata;
   
   assert(scip != NULL);
   assert(cons != NULL);
   assert(violation != NULL);
   
   SCIP_CALL( computeViolation(scip, cons, sol, TRUE) ); /* we assume that scaling was left on */

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   *violation = MAX(consdata->lhsviol, consdata->rhsviol);
   
   return SCIP_OKAY;
}

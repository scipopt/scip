/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   prop_distance.c
 * @ingroup DEFPLUGINS_PROP
 * @brief  propagator for Euclidean distance constraints
 * @author Christopher Hojny
 *
 * This propagator tightens variable domains based on Euclidean distance constraints \f$\| y-z \|_2^2 \geq \delta^2\f$,
 * where both \f$y\f$ and \f$z\f$ are \f$d\f$-dimensional variable vectors and \f$\delta \geq 0\f$ is a variable or
 * fixed number. These constraints are referred to as minimum distance constraints (mindconss) and, if one of the
 * variable vectors is in fact a fixed vector, as minimum distance w.r.t. fixed point (minfd). Moreover, reductions
 * based on pairs of mindconss and minfdconss are derived, which are called minpd and minfpd, respectively.
 *
 * The methods implemented in this propagator are described in the preprint
 *
 * A computational comparison of handling distance constraints in MINLP@n
 * Christopher Hojny and Leo Liberti,@n
 * https://arxiv.org/abs/2605.02305 (2026)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <math.h>

#include "scip/cons_nonlinear.h"
#include "scip/prop_distance.h"

/* fundamental propagator properties */
#define PROP_NAME              "distance"
#define PROP_DESC              "propagator for Euclidean distance constraints"
#define PROP_PRIORITY                 0 /**< propagator priority */
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP /**< propagation timing mask */

#define DEFAULT_MAXNDCONSSFORPAIRS   20 /**< maximum number of distance constraints for which propagation of pairs is applied */
#define DEFAULT_MULTBISECTION       0.1 /**< multiplier used in bisection to split interval [l,u] at point l + multiplier*(u - l) */
#define DEFAULT_GAPBISECTION      0.001 /**< gap between lower and upper bound at which bisection is stopped */
#define DEFAULT_MAXNITERBISECTION     3 /**< maximum number of iterations used in bisection (0: unbounded) */

#define EVENTHDLR_NAME         "propdistance"
#define EVENTHDLR_DESC         "bound change eventhdlr for distance constraints"

#define EPS                       0.0001
#define NPROPS                        4

/*
 * structures and typedefs
 */

/* enum of different types of distance constraints */
enum DCONS_Type
{
   DCONS_TYPE_MIN          = 0,              /**< minimum distance constraints */
   DCONS_TYPE_MINFIXED     = 1,              /**< minimum distance constraints w.r.t. fixed point */
   DCONS_TYPE_MINPAIR      = 2,              /**< pairs of minimum distance constraints */
   DCONS_TYPE_MINFIXEDPAIR = 3               /**< pairs of minimum distance constraints w.r.t. fixed point */
};
typedef enum DCONS_Type DCONS_TYPE;


/** event data */
struct SCIP_EventData
{
   DCONS_TYPE            type;               /**< type of the distance constraint */
   int                   considx;            /**< index of distance constraint in respective group */
};

/** execution method of distance sub-propagator
 *
 *  Searches for domain propagations. The method is called in the node processing loop.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - propdata        : data of (global) distance propagator
 *  - didrun          : pointer to store whether sub-propagator has been called
 *  - infeasible      : pointer to store whether sub-propagator detected cutoff
 *  - nred            : pointer to store number of domain reductions found by sub-propagator
 */
#define SCIP_DECL_DISTCONSPROP(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROPDATA* propdata, SCIP_Bool* didrun, SCIP_Bool* infeasible, int* nred)

/** typedef for sub-propagators of distance constraints */
typedef struct DistconsProp
{
   SCIP_DECL_DISTCONSPROP ((*prop));         /**< propagation method */
} DISTCONSPROP;

/* container for distance constraints of type
 * \f[
 *    \| x - y \|_2^2 \geq r \cdot \delta + c
 * \f]
 * or
 * \f[
 *    \| x - y \|_2^2 \geq r \cdot \delta^2 + c
 * \f]
 * where \f$x\f$ and \f$y\f$ are lists of variables, \f$r, c\f$ are scalars, and \f$\delta\f$ is a variable or scalar
 */
typedef struct DISTCONS_Min
{
   SCIP_VAR**            vars1;              /**< array of first variables in distance constraint */
   SCIP_VAR**            vars2;              /**< array of second variables in distance constraint */
   int                   lenvars;            /**< length of variable arrays */
   SCIP_VAR*             distvar;            /**< variable modeling the distance \f$\delta\f$ (or NULL) */
   SCIP_Real             distval;            /**< fixed distance \f$c\f$ (if distvar is NULL) */
   SCIP_Real             scalar;             /**< multiplicator \f$r\f$ of distance */
   SCIP_Bool             hasvardist;         /**< whether the constraint has a variable distance */
   SCIP_Bool             issquared;          /**< whether distance constraint is squared */
   SCIP_Longint          lastpropnode;       /**< ID of latest node at which constraint has been propagated */
} DISTCONS_MIN;

/* container for distance constraints of type
 * \f[
 *    \| x - p \|_2^2 \geq r \cdot \delta + c
 * \f]
 * where \f$x\f$ is a list of variables, \f$p\f$ is a fixed point, \f$r, c\f$ are scalars, and \f$\delta\f$ is a variable
 */
typedef struct DISTCONS_MinFixed
{
   SCIP_VAR**            vars;               /**< array of variables in distance */
   SCIP_Real*            point;              /**< array of point in distance */
   int                   lenvars;            /**< length of variable arrays */
   SCIP_VAR*             distvar;            /**< variable modeling the distance \f$\delta\f$ (or NULL) */
   SCIP_Real             distval;            /**< fixed distance \f$c\f$ (if distvar is NULL) */
   SCIP_Real             scalar;             /**< multiplicator \f$r\f$ of distance */
   SCIP_Bool             hasvardist;         /**< whether the constraint has a variable distance */
   SCIP_Longint          lastpropnode;       /**< ID of latest node at which constraint has been propagated */
} DISTCONS_MINFIXED;

/** container for pairs of overlapping minimum distance constraints */
typedef struct DISTCONS_MinPair
{
   int                   dimension;          /**< dimension of distance constraints' ambient space */
   SCIP_VAR*             commonvars[3];      /**< common variables of pair of distance constraints */
   SCIP_VAR*             vars1[3];           /**< non-common variables of first distance constraint */
   SCIP_VAR*             vars2[3];           /**< non-common variables of second distance constraint */
   SCIP_VAR*             distvars[2];        /**< distance variables of distance constraints */
   SCIP_Real             scalars[2];         /**< multipliers of distance variables */
   SCIP_Real             constants[2];       /**< constants of distance constraints */
   SCIP_Bool             issquared[2];       /**< whether distance variables are squared */
   DISTCONS_MIN*         dconss[2];          /**< distance constraints of this pair */
} DISTCONS_MINPAIR;

/** container for pairs of overlapping minimum distance constraints w.r.t. a fixed point */
typedef struct DISTCONS_MinFPair
{
   int                   dimension;          /**< dimension of distance constraints' ambient space */
   SCIP_VAR*             vars[3];            /**< common variables of pair of distance constraints */
   SCIP_Real             point1[3];          /**< fixed point of first distance constraint */
   SCIP_Real             point2[3];          /**< fixed point of first distance constraint */
   SCIP_VAR*             distvars[2];        /**< distance variables of distance constraints */
   SCIP_Real             scalars[2];         /**< multipliers of distance variables */
   SCIP_Real             constants[2];       /**< constants of distance constraints */
   DISTCONS_MINFIXED*    dconss[2];          /**< distance constraints of this pair */
} DISTCONS_MINFPAIR;

/** propagator data */
struct SCIP_PropData
{
   /* information about minimum distance (mind) and minimum distance w.r.t. point (minfd) constraints */
   DISTCONS_MIN**        mindconss;          /**< minimum distance constraints */
   int                   nmindconss;         /**< number of minimum distance constraints */
   int                   lenmindconss;       /**< length of mindistconss */
   DISTCONS_MINFIXED**   minfdconss;         /**< minimum distance constraints w.r.t. a fixed point */
   int                   nminfdconss;        /**< number of minimum distance constraints w.r.t. a fixed point */
   int                   lenminfdconss;      /**< length of minfixeddistconss */
   SCIP_Bool             detectedconss;      /**< whether constraints have already been tried to be detected */

   /* information about pairs of such constraints (minpd and minfpd) */
   int                   maxndconssforpairs; /**< maximum number of distance constraints for which propagation of
                                              *   pairs is applied */
   DISTCONS_MINPAIR**    minpdconss;         /**< minpd constraints */
   int                   nminpdconss;        /**< number of minpd constraints */
   int                   lenminpdconss;      /**< length of minpdconss */
   DISTCONS_MINFPAIR**   minfpdconss;        /**< minfpd constraints */
   int                   nminfpdconss;       /**< number of minfpd constraints */
   int                   lenminfpdconss;     /**< length of minfpdconss */

   /* propagation methods */
   DISTCONSPROP*         props[NPROPS];      /**< propagation methods */
   int                   nprops;             /**< number of propagation methods */

   /* event handler related data to more efficiently propagate bounds */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler to handle bound change events */
   SCIP_EVENTDATA*       eventdatamind;      /**< event data for mind constraints */
   SCIP_EVENTDATA*       eventdataminfd;     /**< event data for minfd constraints */
   SCIP_EVENTDATA*       eventdataminpd;     /**< event data for minpd constraints */
   SCIP_EVENTDATA*       eventdataminfpd;    /**< event data for minfpd constraints */
   SCIP_Bool*            dopropmind;         /**< array storing which mind constraint shall be propagated */
   SCIP_Bool*            dopropminfd;        /**< array storing which minfd constraints shall be propagated */
   SCIP_Bool*            dopropminpd;        /**< array storing which minpd constraints shall be propagated */
   SCIP_Bool*            dopropminfpd;       /**< array storing which minfpd constraints shall be propagated */

   /* general limits */
   SCIP_Real             multbisection;      /**< multiplier used in bisection to split interval [l,u] at point l + multiplier*(u - l) */
   SCIP_Real             gapbisection;       /**< gap between lower and upper bound at which bisection is stopped */
   int                   maxniterbisection;  /**< maximum number of iterations used in bisection (0: unbounded) */
};

/* struct to store information about box domains (up to ambient dimension of 3) */
typedef struct Box
{
   int                   dim;                /**< dimension of ambient space of box */
   SCIP_Real             vertices[24];       /**< vertices of box stored as consecutive list */
   int                   nvertices;          /**< number of vertices */
   SCIP_Real             lbs[3];             /**< lower bounds on box coordinates */
   SCIP_Real             ubs[3];             /**< upper bounds on box coordinates */
} BOX;

/*
 * methods for event handler
 */

/** catches variable events of distance constraints */
static
SCIP_RETCODE catchVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to handle bound change events */
   DISTCONS_MIN**        mindconss,          /**< mind conss */
   int                   nmindconss,         /**< number of mind conss */
   SCIP_Bool**           dopropmind,         /**< array storing which mind conss shall be propagated */
   SCIP_EVENTDATA**      eventdatamind,      /**< event data for mind conss */
   DISTCONS_MINFIXED**   minfdconss,         /**< minfd conss */
   int                   nminfdconss,        /**< number of minfd conss */
   SCIP_Bool**           dopropminfd,        /**< array storing which minfd conss shall be propagated */
   SCIP_EVENTDATA**      eventdataminfd,     /**< event data for minfd conss */
   DISTCONS_MINPAIR**    minpdconss,         /**< minpd conss */
   int                   nminpdconss,        /**< number of minpd conss */
   SCIP_Bool**           dopropminpd,        /**< array storing which minpd conss shall be propagated */
   SCIP_EVENTDATA**      eventdataminpd,     /**< event data for minpd conss */
   DISTCONS_MINFPAIR**   minfpdconss,        /**< minfpd conss */
   int                   nminfpdconss,       /**< number of minfpd conss */
   SCIP_Bool**           dopropminfpd,       /**< array storing which minfpd conss shall be propagated */
   SCIP_EVENTDATA**      eventdataminfpd     /**< event data for minfpd conss */
   )
{
   SCIP_EVENTDATA* eventdata;
   SCIP_EVENTTYPE eventtype;
   int dim;
   int c;
   int d;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(dopropmind != NULL);
   assert(dopropminfd != NULL);
   assert(dopropminpd != NULL);
   assert(dopropminfpd != NULL);
   assert(eventdatamind != NULL);
   assert(eventdataminfd != NULL);
   assert(eventdataminpd != NULL);
   assert(eventdataminfpd != NULL);

   eventtype = SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_UBTIGHTENED;

   if( nmindconss > 0 )
   {
      assert(mindconss != NULL);

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, dopropmind, nmindconss) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, eventdatamind, nmindconss) );
      for( c = 0; c < nmindconss; ++c )
      {
         assert(mindconss[c] != NULL);

         eventdata = &(*eventdatamind)[c];
         eventdata->type = DCONS_TYPE_MIN;
         eventdata->considx = c;

         dim = mindconss[c]->lenvars;
         for( d = 0; d < dim; ++d )
         {
            SCIP_CALL( SCIPcatchVarEvent(scip, mindconss[c]->vars1[d], eventtype, eventhdlr, eventdata, NULL) );
            SCIP_CALL( SCIPcatchVarEvent(scip, mindconss[c]->vars2[d], eventtype, eventhdlr, eventdata, NULL) );
         }
         if( mindconss[c]->distvar != NULL )
         {
            SCIP_CALL( SCIPcatchVarEvent(scip, mindconss[c]->distvar, SCIP_EVENTTYPE_LBTIGHTENED,
                  eventhdlr, eventdata, NULL) );
         }

         (*dopropmind)[c] = FALSE;
      }
   }
   if( nminfdconss > 0 )
   {
      assert(minfdconss != NULL);

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, dopropminfd, nminfdconss) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, eventdataminfd, nminfdconss) );
      for( c = 0; c < nminfdconss; ++c )
      {
         assert(minfdconss[c] != NULL);

         eventdata = &(*eventdataminfd)[c];
         eventdata->type = DCONS_TYPE_MINFIXED;
         eventdata->considx = c;

         dim = minfdconss[c]->lenvars;
         for( d = 0; d < dim; ++d )
         {
            SCIP_CALL( SCIPcatchVarEvent(scip, minfdconss[c]->vars[d], eventtype, eventhdlr, eventdata, NULL) );
         }
         if( minfdconss[c]->distvar != NULL )
         {
            SCIP_CALL( SCIPcatchVarEvent(scip, minfdconss[c]->distvar, SCIP_EVENTTYPE_LBTIGHTENED,
                  eventhdlr, eventdata, NULL) );
         }

         (*dopropminfd)[c] = FALSE;
      }
   }
   if( nminpdconss > 0 )
   {
      assert(minpdconss != NULL);

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, dopropminpd, nminpdconss) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, eventdataminpd, nminpdconss) );
      for( c = 0; c < nminpdconss; ++c )
      {
         assert(minpdconss[c] != NULL);

         eventdata = &(*eventdataminpd)[c];
         eventdata->type = DCONS_TYPE_MINPAIR;
         eventdata->considx = c;

         dim = minpdconss[c]->dimension;
         for( d = 0; d < dim; ++d )
         {
            SCIP_CALL( SCIPcatchVarEvent(scip, minpdconss[c]->commonvars[d], eventtype, eventhdlr, eventdata, NULL) );
            SCIP_CALL( SCIPcatchVarEvent(scip, minpdconss[c]->vars1[d], eventtype, eventhdlr, eventdata, NULL) );
            SCIP_CALL( SCIPcatchVarEvent(scip, minpdconss[c]->vars2[d], eventtype, eventhdlr, eventdata, NULL) );
         }
         if( minpdconss[c]->distvars[0] != NULL )
         {
            SCIP_CALL( SCIPcatchVarEvent(scip, minpdconss[c]->distvars[0], SCIP_EVENTTYPE_LBTIGHTENED,
                  eventhdlr, eventdata, NULL) );
         }
         if( minpdconss[c]->distvars[1] != NULL )
         {
            SCIP_CALL( SCIPcatchVarEvent(scip, minpdconss[c]->distvars[1], SCIP_EVENTTYPE_LBTIGHTENED,
                  eventhdlr, eventdata, NULL) );
         }
         (*dopropminpd)[c] = FALSE;
      }
   }
   if( nminfpdconss > 0 )
   {
      assert(minfpdconss != NULL);

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, dopropminfpd, nminfpdconss) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, eventdataminfpd, nminfpdconss) );
      for( c = 0; c < nminfpdconss; ++c )
      {
         assert(minfpdconss[c] != NULL);

         eventdata = &(*eventdataminfpd)[c];
         eventdata->type = DCONS_TYPE_MINFIXEDPAIR;
         eventdata->considx = c;

         dim = minfpdconss[c]->dimension;
         for( d = 0; d < dim; ++d )
         {
            SCIP_CALL( SCIPcatchVarEvent(scip, minfpdconss[c]->vars[d], eventtype, eventhdlr, eventdata, NULL) );
         }
         if( minfpdconss[c]->distvars[0] != NULL )
         {
            SCIP_CALL( SCIPcatchVarEvent(scip, minfpdconss[c]->distvars[0], SCIP_EVENTTYPE_LBTIGHTENED,
                  eventhdlr, eventdata, NULL) );
         }
         if( minfpdconss[c]->distvars[1] != NULL )
         {
            SCIP_CALL( SCIPcatchVarEvent(scip, minfpdconss[c]->distvars[1], SCIP_EVENTTYPE_LBTIGHTENED,
                  eventhdlr, eventdata, NULL) );
         }
         (*dopropminfpd)[c] = FALSE;
      }
   }

   return SCIP_OKAY;
}


/** drops variable events of distance constraints */
static
SCIP_RETCODE dropVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to handle bound change events */
   DISTCONS_MIN**        mindconss,          /**< mind conss */
   int                   nmindconss,         /**< number of mind conss w.r.t. a fixed point */
   SCIP_EVENTDATA*       eventdatamind,      /**< event data for mind conss */
   DISTCONS_MINFIXED**   minfdconss,         /**< minfd conss */
   int                   nminfdconss,        /**< number of minfd conss */
   SCIP_EVENTDATA*       eventdataminfd,     /**< event data for minfd conss */
   DISTCONS_MINPAIR**    minpdconss,         /**< minpd conss */
   int                   nminpdconss,        /**< number of minpd conss */
   SCIP_EVENTDATA*       eventdataminpd,     /**< event data for minpd conss */
   DISTCONS_MINFPAIR**   minfpdconss,        /**< minfpd conss */
   int                   nminfpdconss,       /**< number of minfpd conss */
   SCIP_EVENTDATA*       eventdataminfpd     /**< event data for minfpd conss */
   )
{
   SCIP_EVENTTYPE eventtype;
   int dim;
   int c;
   int d;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventtype = SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_UBTIGHTENED;

   if( nmindconss > 0 )
   {
      assert(mindconss != NULL);

      for( c = 0; c < nmindconss; ++c )
      {
         assert(mindconss[c] != NULL);

         dim = mindconss[c]->lenvars;
         for( d = 0; d < dim; ++d )
         {
            SCIP_CALL( SCIPdropVarEvent(scip, mindconss[c]->vars1[d], eventtype, eventhdlr, &eventdatamind[c], -1) );
            SCIP_CALL( SCIPdropVarEvent(scip, mindconss[c]->vars2[d], eventtype, eventhdlr, &eventdatamind[c], -1) );
         }
         if( mindconss[c]->distvar != NULL )
         {
            SCIP_CALL( SCIPdropVarEvent(scip, mindconss[c]->distvar, SCIP_EVENTTYPE_LBTIGHTENED,
                  eventhdlr, &eventdatamind[c], -1) );
         }
      }
   }
   if( nminfdconss > 0 )
   {
      assert(minfdconss != NULL);

      for( c = 0; c < nminfdconss; ++c )
      {
         assert(minfdconss[c] != NULL);

         dim = minfdconss[c]->lenvars;
         for( d = 0; d < dim; ++d )
         {
            SCIP_CALL( SCIPdropVarEvent(scip, minfdconss[c]->vars[d], eventtype, eventhdlr, &eventdataminfd[c], -1) );
         }
         if( minfdconss[c]->distvar != NULL )
         {
            SCIP_CALL( SCIPdropVarEvent(scip, minfdconss[c]->distvar, SCIP_EVENTTYPE_LBTIGHTENED,
                  eventhdlr, &eventdataminfd[c], -1) );
         }
      }
   }
   if( nminpdconss > 0 )
   {
      assert(minpdconss != NULL);

      for( c = 0; c < nminpdconss; ++c )
      {
         assert(minpdconss[c] != NULL);

         dim = minpdconss[c]->dimension;
         for( d = 0; d < dim; ++d )
         {
            SCIP_CALL( SCIPdropVarEvent(scip, minpdconss[c]->commonvars[d], eventtype, eventhdlr,
                  &eventdataminpd[c], -1) );
            SCIP_CALL( SCIPdropVarEvent(scip, minpdconss[c]->vars1[d], eventtype, eventhdlr, &eventdataminpd[c], -1) );
            SCIP_CALL( SCIPdropVarEvent(scip, minpdconss[c]->vars2[d], eventtype, eventhdlr, &eventdataminpd[c], -1) );
         }
         if( minpdconss[c]->distvars[0] != NULL )
         {
            SCIP_CALL( SCIPdropVarEvent(scip, minpdconss[c]->distvars[0], SCIP_EVENTTYPE_LBTIGHTENED,
                  eventhdlr, &eventdataminpd[c], -1) );
         }
         if( minpdconss[c]->distvars[1] != NULL )
         {
            SCIP_CALL( SCIPdropVarEvent(scip, minpdconss[c]->distvars[1], SCIP_EVENTTYPE_LBTIGHTENED,
                  eventhdlr, &eventdataminpd[c], -1) );
         }
      }
   }
   if( nminfpdconss > 0 )
   {
      assert(minfpdconss != NULL);

      for( c = 0; c < nminfpdconss; ++c )
      {
         assert(minfpdconss[c] != NULL);

         dim = minfpdconss[c]->dimension;
         for( d = 0; d < dim; ++d )
         {
            SCIP_CALL( SCIPdropVarEvent(scip, minfpdconss[c]->vars[d], eventtype, eventhdlr, &eventdataminfpd[c], -1) );
         }
         if( minfpdconss[c]->distvars[0] != NULL )
         {
            SCIP_CALL( SCIPdropVarEvent(scip, minfpdconss[c]->distvars[0], SCIP_EVENTTYPE_LBTIGHTENED,
                  eventhdlr, &eventdataminfpd[c], -1) );
         }
         if( minfpdconss[c]->distvars[1] != NULL )
         {
            SCIP_CALL( SCIPdropVarEvent(scip, minfpdconss[c]->distvars[1], SCIP_EVENTTYPE_LBTIGHTENED,
                  eventhdlr, &eventdataminfpd[c], -1) );
         }
      }
   }

   return SCIP_OKAY;
}


/** execution method of bound change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecDistcons)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);

   propdata = (SCIP_PROPDATA*)SCIPeventhdlrGetData(eventhdlr);

   switch( eventdata->type )
   {
   case DCONS_TYPE_MIN:
      propdata->dopropmind[eventdata->considx] = TRUE;
      break;
   case DCONS_TYPE_MINFIXED:
      propdata->dopropminfd[eventdata->considx] = TRUE;
      break;
   case DCONS_TYPE_MINPAIR:
      propdata->dopropminpd[eventdata->considx] = TRUE;
      break;
   default:
      assert(eventdata->type == DCONS_TYPE_MINFIXEDPAIR);
      propdata->dopropminfpd[eventdata->considx] = TRUE;
   }

   return SCIP_OKAY;
}

/*
 * methods for distance constraints
 */

/** tries to extract a minfd constraint from expression based on squared differences */
static
SCIP_RETCODE tryExtractDistconsMinfd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression modeling nonlinear constraint */
   SCIP_Real             lhs,                /**< left-hand side of constraint */
   SCIP_Real             rhs,                /**< right-hand side of constraint */
   DISTCONS_MINFIXED***  distconss,          /**< pointer to store minfd conss */
   int*                  ndistconss,         /**< pointer to store number of minfd */
   int*                  lendistconss        /**< pointer to store length of minfd conss */
   )
{
   SCIP_Real constant = 0.0;
   SCIP_EXPR** children;
   SCIP_VAR** powvars;
   SCIP_VAR** linvars;
   SCIP_VAR* distvar = NULL;
   int* pow2lin = NULL;
   SCIP_Real* lincoefs;
   SCIP_Real distcoef = 0.0;
   SCIP_Real powcoef;
   SCIP_Bool isleq;
   SCIP_Bool success;
   int nchildren;
   int npowvars = 0;
   int nlinvars = 0;
   int start;
   int i;
   int j;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(distconss != NULL);
   assert(ndistconss != NULL);
   assert(lendistconss != NULL);

   success = FALSE;

   /* currently, we only handle inequalities with one infinite bound */
   if( !SCIPisInfinity(scip, rhs) && !SCIPisInfinity(scip, -lhs) )
      return SCIP_OKAY;

   /* the root expression needs to be a sum expression */
   if( !SCIPisExprSum(scip, expr) )
      return SCIP_OKAY;

   /* determine type of constraint and its expected coefficients */
   assert(!SCIPisInfinity(scip, lhs) ||!SCIPisInfinity(scip, rhs));
   isleq = SCIPisInfinity(scip, -lhs);
   powcoef = isleq ? -1.0 : 1.0;

   constant = SCIPgetConstantExprSum(expr);
   nchildren = SCIPexprGetNChildren(expr);
   children = SCIPexprGetChildren(expr);

   /* extract information about variables, powers, and constants */
   SCIP_CALL( SCIPallocBufferArray(scip, &powvars, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvars, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, nchildren) );
   for( i = 0; i < nchildren; ++i )
   {
      SCIP_EXPR* child;

      child = children[i];

      if( SCIPisExprVar(scip, child) )
      {
         linvars[nlinvars] = SCIPgetVarExprVar(child);
         lincoefs[nlinvars++] = SCIPgetCoefsExprSum(expr)[i];
      }
      else if( SCIPisExprPower(scip, child) )
      {
         SCIP_Real pow;

         pow = SCIPgetExponentExprPow(child);
         if( !SCIPisEQ(scip, pow, 2.0) )
            goto FREEMEMORY;
         if( !SCIPisExprVar(scip, SCIPexprGetChildren(child)[0]) )
            goto FREEMEMORY;

         if( !SCIPisEQ(scip, SCIPgetCoefsExprSum(expr)[i], powcoef) )
            goto FREEMEMORY;
         else
            powvars[npowvars++] = SCIPgetVarExprVar(SCIPexprGetChildren(child)[0]);
      }
      else
      {
         assert(!SCIPisExprValue(scip, child));
         goto FREEMEMORY;
      }
   }

   /* in a distance constraint, every linear variable has a power variable (except for at most one variable) */
   SCIP_CALL( SCIPallocBufferArray(scip, &pow2lin, npowvars) );
   for( i = 0; i < npowvars; ++i )
      pow2lin[i] = -1;         /* map power variables -> linear variables */
   success = TRUE;             /* needed in case there is no linear variable */
   start = 0;
   for( i = 0; i < nlinvars && success; ++i )
   {
      success = FALSE;
      for( j = start; j < npowvars; ++j )
      {
         if( pow2lin[j] != -1 )
            continue;

         if( linvars[i] == powvars[j] )
         {
            success = TRUE;
            pow2lin[j] = i;
            if( j == start )
               ++start;
         }
      }

      /* one unmatched variable can be turned into a distance variable */
      if( !success && distvar == NULL )
      {
         distvar = linvars[i];
         distcoef = lincoefs[i];
         success = TRUE;
      }
   }

   if( !success )
      goto FREEMEMORY;

   /* store constraint */
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, distconss, lendistconss, *ndistconss + 1) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &(*distconss)[*ndistconss]) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*distconss)[*ndistconss]->vars, npowvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*distconss)[*ndistconss]->point, npowvars) );
   for( j = 0; j < npowvars; ++j )
   {
      /* find scalar p(oint) such that (var - p)^2 = var^2 - 2*p*var + p^2 */
      (*distconss)[*ndistconss]->vars[j] = powvars[j];
      if( pow2lin[j] == -1 )
         (*distconss)[*ndistconss]->point[j] = 0.0;
      else
      {
         SCIP_Real coord;

         coord = lincoefs[pow2lin[j]] / 2;

         /* store point (normalize based on inequality type) */
         (*distconss)[*ndistconss]->point[j] = isleq ? coord : -coord;

         /* remove p^2 from the constant */
         if( isleq )
            constant += coord * coord;
         else
            constant -= coord * coord;
      }
   }

   /* adapt constant based on type of inequality */
   if( isleq )
      constant -= rhs;
   else
      constant = lhs - constant;

   if( distvar != NULL )
   {
      (*distconss)[*ndistconss]->distvar = distvar;
      (*distconss)[*ndistconss]->hasvardist = TRUE;
      (*distconss)[*ndistconss]->scalar = isleq ? distcoef : -distcoef;
      (*distconss)[*ndistconss]->distval = constant;
   }
   else
   {
      (*distconss)[*ndistconss]->distvar = NULL;
      (*distconss)[*ndistconss]->hasvardist = FALSE;
      (*distconss)[*ndistconss]->scalar = 1.0;
      (*distconss)[*ndistconss]->distval = constant;
   }
   (*distconss)[*ndistconss]->lenvars = npowvars;
   (*distconss)[*ndistconss]->lastpropnode = -1;
   (*ndistconss)++;

 FREEMEMORY:
   SCIPfreeBufferArrayNull(scip, &pow2lin);
   SCIPfreeBufferArray(scip, &lincoefs);
   SCIPfreeBufferArray(scip, &linvars);
   SCIPfreeBufferArray(scip, &powvars);

   return SCIP_OKAY;
}

/** tries to extract a mind constraint from a nonlinear constraint */
static
SCIP_RETCODE tryExtractDistconsMind(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   DISTCONS_MIN***       distconss,          /**< pointer to store minimum distance constraints */
   int*                  ndistconss,         /**< pointer to store number of minimum distance constraints */
   int*                  lendistconss        /**< pointer to store length of mindistconss */
   )
{
   SCIP_EXPR** summands;
   SCIP_VAR** powvars;
   SCIP_VAR** prodvars;
   SCIP_EXPR* expr;
   SCIP_Real* sumcoefs;
   SCIP_Bool* used;
   SCIP_Real prodcoef;
   SCIP_Real powcoef;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int ndistterms;
   int nsummands;
   int start;
   int c;
   SCIP_Bool success = FALSE;
   SCIP_VAR* boundingvar = NULL;
   SCIP_Bool issquared = FALSE;
   int npowexprs = 0;
   int nprodexprs = 0;
   SCIP_Real boundingcoef = 0.0;
   SCIP_Bool isleq;
   SCIP_Real constant;

   assert(scip != NULL);

   expr = SCIPgetExprNonlinear(cons);
   lhs = SCIPgetLhsNonlinear(cons);
   rhs = SCIPgetRhsNonlinear(cons);
   assert(expr != NULL);
   assert(!SCIPisInfinity(scip, rhs) || !SCIPisInfinity(scip, -lhs));

   /* currently, we only handle inequalities with one infinite bound */
   if( !SCIPisInfinity(scip, rhs) && !SCIPisInfinity(scip, -lhs) )
      return SCIP_OKAY;

   /*
    * We scan for expressions of type sum(x_{i,d}^2 - 2x_{i,d}x_{j,d} + x_{j,d}^2, i,j,d) + rest,
    * where rest is const*y or const*y^2.
    */

   /* the root expression needs to be a sum */
   if( !SCIPisExprSum(scip, expr) )
      return SCIP_OKAY;

   isleq = SCIPisInfinity(scip, rhs) ? FALSE : TRUE;
   constant = isleq ? -rhs : -lhs; /* interpret the rhs/lhs as part of the expression */
   prodcoef = isleq ? 2.0 : -2.0;
   powcoef = isleq ? -1.0 : 1.0;

   nsummands = SCIPexprGetNChildren(expr);
   summands = SCIPexprGetChildren(expr);
   sumcoefs = SCIPgetCoefsExprSum(expr);

   constant += SCIPgetConstantExprSum(expr);

   ndistterms = 2 * nsummands / 3 + 1;
   SCIP_CALL( SCIPallocBufferArray(scip, &powvars, ndistterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &prodvars, ndistterms) );

   for( c = 0; c < nsummands; ++c )
   {
      SCIP_EXPR* child;

      child = summands[c];

      if( SCIPisExprVar(scip, child) )
      {
         /* possibly store the comparison variable, there is at most one of them */
         if( boundingvar == NULL )
         {
            boundingvar = SCIPgetVarExprVar(child);
            boundingcoef = sumcoefs[c];
         }
         else
            break;
      }
      else if( SCIPisExprPower(scip, child) )
      {
         SCIP_EXPR* powchild;

         powchild = SCIPexprGetChildren(child)[0];
         assert(powchild != NULL);

         /* if there are too many power expressions, it canot be a distance constraint */
         if( npowexprs >= ndistterms )
            break;
         if( !SCIPisExprVar(scip, powchild) )
            break;
         if( !SCIPisEQ(scip, SCIPgetExponentExprPow(child), 2.0) )
            break;
         if( !SCIPisEQ(scip, sumcoefs[c], powcoef) )
         {
            if( boundingvar == NULL )
            {
               boundingvar = SCIPgetVarExprVar(powchild);
               boundingcoef = sumcoefs[c];
               issquared = TRUE;
            }
            else
               break;
         }
         else
            powvars[npowexprs++] = SCIPgetVarExprVar(powchild);
      }
      else if( SCIPisExprProduct(scip, child) )
      {
         /* we can only handle products of two variables */
         if( SCIPexprGetNChildren(child) != 2 )
            break;

         /* if there are too many product expressions, it canot be a distance constraint */
         if( nprodexprs >= ndistterms - 1 )
            break;

         assert(SCIPgetCoefExprProduct(child) == 1.0);
         if( !SCIPisEQ(scip, prodcoef, sumcoefs[c]) )
            break;

         if( !SCIPisExprVar(scip, SCIPexprGetChildren(child)[0])
            || !SCIPisExprVar(scip, SCIPexprGetChildren(child)[1]) )
            break;
         prodvars[nprodexprs++] = SCIPgetVarExprVar(SCIPexprGetChildren(child)[0]);
         prodvars[nprodexprs++] = SCIPgetVarExprVar(SCIPexprGetChildren(child)[1]);
      }
      else
      {
         assert(!SCIPisExprValue(scip, child));
         break;
      }
   }

   /* if we terminated early, we had no success */
   if( c < nsummands )
      goto FREEMEMORY;

   /* if there are no product variables, we cannot detect a distance constraint */
   if( nprodexprs <= 0 )
      goto FREEMEMORY;

   /* if the number of power and product expressions do not match, we cannot detect a distance constraint */
   if( nprodexprs != npowexprs )
      goto FREEMEMORY;

   /* if the minimal distance is constraint to be zero, then it is trivially satisfied */
   if( boundingvar == NULL && SCIPisZero(scip, constant) )
      goto FREEMEMORY;

   /* check whether the product and power expressions match */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &used, nprodexprs) );
   start = 0;
   success = TRUE;
   for( c = 0; c < nprodexprs / 2; ++c )
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_Bool found1 = FALSE;
      SCIP_Bool found2 = FALSE;
      int i;

      var1 = prodvars[2 * c];
      var2 = prodvars[2 * c + 1];

      /* check whether the power variables match the product variables, keep track of the first "start"
       * many variables that have been matched (most likely, the constraints are ordered)
       */
      for( i = start; i < npowexprs && (!found1 || !found2); ++i )
      {
         if( used[i] )
            continue;

         if( !found1 && powvars[i] == var1 )
         {
            used[i] = TRUE;
            found1 = TRUE;
            if( i == start )
               ++start;
         }
         else if( !found2 && powvars[i] == var2 )
         {
            used[i] = TRUE;
            found2 = TRUE;
            if( i == start )
               ++start;
         }
      }

      if( !found1 || !found2 )
      {
         success = FALSE;
         break;
      }
   }
   SCIPfreeBufferArray(scip, &used);

   if( !success )
      goto FREEMEMORY;

   /* ensure that memory suffices to store constraint */
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, distconss, lendistconss, *ndistconss + 1) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &(*distconss)[*ndistconss]) );

   (*distconss)[*ndistconss]->lenvars = nprodexprs / 2;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*distconss)[*ndistconss]->vars1, nprodexprs / 2) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*distconss)[*ndistconss]->vars2, nprodexprs / 2) );
   for( c = 0; c < nprodexprs / 2; ++c )
   {
      (*distconss)[*ndistconss]->vars1[c] = prodvars[2*c];
      (*distconss)[*ndistconss]->vars2[c] = prodvars[2*c + 1];
   }
   (*distconss)[*ndistconss]->distvar = boundingvar;
   (*distconss)[*ndistconss]->issquared = issquared;
   (*distconss)[*ndistconss]->scalar = isleq ? boundingcoef : -boundingcoef;
   (*distconss)[*ndistconss]->distval = isleq ? constant : -constant;
   (*distconss)[*ndistconss]->hasvardist = boundingvar != NULL;
   (*distconss)[*ndistconss]->lastpropnode = -1;
   (*ndistconss)++;

 FREEMEMORY:
   SCIPfreeBufferArray(scip, &prodvars);
   SCIPfreeBufferArray(scip, &powvars);

   return SCIP_OKAY;
}

/** detects mind constraints */
static
SCIP_RETCODE detectDistconssMind(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlrnonlinear,  /**< nonlinear constraint handler */
   DISTCONS_MIN***       distconss,          /**< pointer to store minimum distance constraints */
   int*                  ndistconss,         /**< pointer to store number of minimum distance constraints */
   int*                  lendistconss        /**< pointer to store length of mindistconss */
   )
{
   SCIP_CONS** conss;
   int nconss;
   int c;

   assert(scip != NULL);
   assert(conshdlrnonlinear != NULL);

   conss = SCIPconshdlrGetConss(conshdlrnonlinear);
   nconss = SCIPconshdlrGetNConss(conshdlrnonlinear);
   assert(conss != NULL || nconss == 0);

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( tryExtractDistconsMind(scip, conss[c], distconss, ndistconss, lendistconss) );
   }

   return SCIP_OKAY;
}

/** detects minfd constraints */
static
SCIP_RETCODE detectDistconssMinfd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlrnonlinear,  /**< nonlinear constraint handler */
   DISTCONS_MINFIXED***  distconss,          /**< pointer to store minimum distance constraints w.r.t. a fixed point */
   int*                  ndistconss,         /**< pointer to store number of min. dist. conss w.r.t. a fixed point */
   int*                  lendistconss        /**< pointer to store length of distconss */
   )
{
   SCIP_CONS** conss;
   SCIP_EXPR* expr;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int nconss;
   int c;

   assert(scip != NULL);
   assert(conshdlrnonlinear != NULL);

   conss = SCIPconshdlrGetConss(conshdlrnonlinear);
   nconss = SCIPconshdlrGetNConss(conshdlrnonlinear);
   assert(conss != NULL || nconss == 0);

   for( c = 0; c < nconss; ++c )
   {
      expr = SCIPgetExprNonlinear(conss[c]);
      lhs = SCIPgetLhsNonlinear(conss[c]);
      rhs = SCIPgetRhsNonlinear(conss[c]);

      SCIP_CALL( tryExtractDistconsMinfd(scip, expr, lhs, rhs, distconss, ndistconss, lendistconss) );
   }

   return SCIP_OKAY;
}

/** extracts minimum distance constraints from problem data */
static
SCIP_RETCODE getMinDistanceConss(
   SCIP*                 scip,               /**< SCIP data structure */
   DISTCONS_MIN***       mindconss,          /**< pointer to store mind conss */
   int*                  nmindconss,         /**< pointer to store number of mind conss */
   int*                  lenmindconss,       /**< pointer to store length of mindconss */
   DISTCONS_MINFIXED***  minfdconss,         /**< pointer to store minfd conss */
   int*                  nminfdconss,        /**< pointer to store number of minfd conss */
   int*                  lenminfdconss       /**< pointer to store length of minfdconss */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(scip != NULL);

   /* if the nonlinear constraint handler is not present, we cannot detect minimum distance constraints */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   if( conshdlr == NULL )
      return SCIP_OKAY;

   SCIP_CALL( detectDistconssMind(scip, conshdlr, mindconss, nmindconss, lenmindconss) );
   SCIP_CALL( detectDistconssMinfd(scip, conshdlr, minfdconss, nminfdconss, lenminfdconss) );

   return SCIP_OKAY;
}

/** frees a mind constraint */
static
SCIP_RETCODE freeDistconsMind(
   SCIP*                 scip,               /**< SCIP data structure */
   DISTCONS_MIN**        distcons            /**< pointer to mind constraint */
   )
{
   assert(scip != NULL);
   assert(distcons != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &(*distcons)->vars1, (*distcons)->lenvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*distcons)->vars2, (*distcons)->lenvars);

   SCIPfreeBlockMemory(scip, distcons);

   return SCIP_OKAY;
}

/** frees a minfd constraint */
static
SCIP_RETCODE freeDistconsMinfd(
   SCIP*                 scip,               /**< SCIP data structure */
   DISTCONS_MINFIXED**   distcons            /**< pointer to minfd constraint */
   )
{
   assert(scip != NULL);
   assert(distcons != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &(*distcons)->vars, (*distcons)->lenvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*distcons)->point, (*distcons)->lenvars);

   SCIPfreeBlockMemory(scip, distcons);

   return SCIP_OKAY;
}

/*
 * general functions
 */

/** returns maximum distance between two variables */
static
SCIP_Real maxDistVarPair(
   SCIP_VAR*             var1,               /**< first variable in variable pair */
   SCIP_VAR*             var2                /**< second variable in variable pair */
   )
{
   SCIP_Real val1;
   SCIP_Real val2;

   assert(var1 != NULL);
   assert(var2 != NULL);

   val1 = SCIPvarGetUbLocal(var1) - SCIPvarGetLbLocal(var2);
   val2 = SCIPvarGetUbLocal(var2) - SCIPvarGetLbLocal(var1);

   return MAX(val1, val2);
}

/** tries to improve a variable's lower bound */
static
SCIP_RETCODE improveLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose bound shall be improved */
   SCIP_Real             bound,              /**< new bound */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility has been detected */
   int*                  nred                /**< pointer to increment by number of found reductions */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(infeasible != NULL);
   assert(nred != NULL);

   SCIP_Bool tightened;

   assert(scip != NULL);
   assert(var != NULL);
   assert(infeasible != NULL);
   assert(nred != NULL);

   SCIP_CALL( SCIPtightenVarLb(scip, var, bound, FALSE, infeasible, &tightened) );

   if( tightened )
      ++(*nred);

   return SCIP_OKAY;
}

/** tries to improve a variable's upper bound */
static
SCIP_RETCODE improveUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose bound shall be improved */
   SCIP_Real             bound,              /**< new bound */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility has been detected */
   int*                  nred                /**< pointer to increment by number of found reductions */
   )
{
   SCIP_Bool tightened;

   assert(scip != NULL);
   assert(var != NULL);
   assert(infeasible != NULL);
   assert(nred != NULL);

   SCIP_CALL( SCIPtightenVarUb(scip, var, bound, FALSE, infeasible, &tightened) );

   if( tightened )
      ++(*nred);

   return SCIP_OKAY;
}

/** returns the minimum value of expression \f$a \cdot x\f$ or \f$a \cdot x^2\f$ for value of \f$x\f$ in an interval,
 *  where \f$x\f$ corresponds to a distance variable */
static
SCIP_Real getMinimumScalarVariableProduct(
   SCIP_Real             lb,                 /**< lower bound on variable */
   SCIP_Real             ub,                 /**< upper bound on variable */
   SCIP_Real             scalar,             /**< scalar in front of variable */
   SCIP_Bool             issquared           /**< whether we consider a squared variable expression */
   )
{
   if( !issquared )
      return scalar < 0.0 ? scalar * ub : scalar * lb;

   if( scalar < 0.0 )
      return scalar * MIN(lb * lb, ub * ub);

   if( lb > 0.0 )
      return scalar * lb * lb;

   if( ub < 0.0 )
      return scalar * ub * ub;

   return 0.0;
}

/** returns the squared ell2-distance of two points */
static
SCIP_Real getDistL2Sq(
   SCIP_Real*            point1,             /**< first point */
   SCIP_Real*            point2,             /**< second point */
   int                   dim                 /**< length of points */
   )
{
   SCIP_Real res = 0.0;
   int d;

   assert(point1 != NULL);
   assert(point2 != NULL);

   for( d = 0; d < dim; ++d )
      res += SQR(point1[d] - point2[d]);

   return res;
}

/*
 * methods for propagating different types of distance constraints
 */

/** propagates a 1-dimensional component of mind cons using a simple propagation method,
 *  see Proposition 2.3 in the preprint
 *  C. Hojny and L. Liberti, A computational comparison of handling distance constraints in MINLP, 2026. */
static
SCIP_RETCODE propagateMindConsSimple1D(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var1,               /**< first variable in variable pair */
   SCIP_VAR*             var2,               /**< second variable in variable pair */
   SCIP_Real             sqdistbound,        /**< squared distance bound for propagation */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility has been detected */
   int*                  nred                /**< pointer to increment by number of found reductions */
   )
{
   SCIP_Real distbound;
   SCIP_Real lb1;
   SCIP_Real lb2;
   SCIP_Real ub1;
   SCIP_Real ub2;

   assert(scip != NULL);
   assert(var1 != NULL);
   assert(var2 != NULL);
   assert(infeasible != NULL);
   assert(nred != NULL);

   if( sqdistbound < 0.0 )
      return SCIP_OKAY;

   lb1 = SCIPvarGetLbLocal(var1);
   lb2 = SCIPvarGetLbLocal(var2);
   ub1 = SCIPvarGetUbLocal(var1);
   ub2 = SCIPvarGetUbLocal(var2);

   distbound = sqrt(sqdistbound);

   if( SCIPisLE(scip, ub1, lb2) )
   {
      /* var1 is to the left of var2 */
      SCIP_CALL( improveUb(scip, var1, ub2 - distbound, infeasible, nred) );
      if( !*infeasible )
      {
         SCIP_CALL( improveLb(scip, var2, lb1 + distbound, infeasible, nred) );
      }
   }
   else if( SCIPisGE(scip, lb1, ub2) )
   {
      /* var2 is to the left of var1 */
      SCIP_CALL( improveUb(scip, var2, ub1 - distbound, infeasible, nred) );
      if( !*infeasible )
      {
         SCIP_CALL( improveLb(scip, var1, lb2 + distbound, infeasible, nred) );
      }
   }
   else if( SCIPisGE(scip, ub1, ub2) && SCIPisLE(scip, lb1, lb2) )
   {
      /* var1's domain contains var2's domain */
      if( SCIPisLT(scip, ub2 - lb1, distbound) && SCIPisLT(scip, ub1 - lb2, distbound) )
         *infeasible = TRUE;
   }
   else if( SCIPisGE(scip, ub2, ub1) && SCIPisLE(scip, lb2, lb1) )
   {
      /* var2's domain contains var1's domain */
      if( SCIPisLT(scip, ub1 - lb2, distbound) && SCIPisLT(scip, ub2 - lb1, distbound) )
         *infeasible = TRUE;
   }
   else if( SCIPisLT(scip, lb1, lb2) )
   {
      /* both domains overlap and var1 is further to the left */
      if( SCIPisLT(scip, ub1 - lb2, distbound) && SCIPisLT(scip, ub2 - ub1, distbound) )
      {
         SCIP_CALL( improveUb(scip, var1, ub2 - distbound, infeasible, nred) );
      }
      if( !*infeasible && SCIPisLT(scip, lb2 - lb1, distbound) && SCIPisLT(scip, ub1 - lb2, distbound) )
      {
         SCIP_CALL( improveLb(scip, var2, lb1 + distbound, infeasible, nred) );
      }
   }
   else if( SCIPisLT(scip, lb2, lb1) )
   {
      /* both domains overlap and var2 is further to the left */
      if( SCIPisLT(scip, ub2 - lb1, distbound) && SCIPisLT(scip, ub1 - ub2, distbound) )
      {
         SCIP_CALL( improveUb(scip, var2, ub1 - distbound, infeasible, nred) );
      }
      if( !*infeasible && SCIPisLT(scip, lb1 - lb2, distbound) && SCIPisLT(scip, ub2 - lb1, distbound) )
      {
         SCIP_CALL( improveLb(scip, var1, lb2 + distbound, infeasible, nred) );
      }
   }

   return SCIP_OKAY;
}


/** propagates a mind cons using a simple propagation method
 *  see Section 2.1.2 in the preprint
 *  C. Hojny and L. Liberti, A computational comparison of handling distance constraints in MINLP, 2026. */
static
SCIP_RETCODE propagateMindConsSimple(
   SCIP*                 scip,               /**< SCIP data structure */
   DISTCONS_MIN*         mindcons,           /**< mind cons */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility is detected */
   int*                  nred                /**< pointer to store number of found reductions */
   )
{
   SCIP_Real totalsquareddist = 0.0;
   SCIP_Real distbound = 0.0;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   SCIP_Real dist;
   int nvars;
   int i;

   assert(mindcons != NULL);
   assert(infeasible != NULL);
   assert(nred != NULL);

   *nred = 0;
   *infeasible = FALSE;

   /* get data */
   vars1 = mindcons->vars1;
   vars2 = mindcons->vars2;
   nvars = mindcons->lenvars;
   if( mindcons->hasvardist )
   {
      SCIP_Real lb;
      SCIP_Real ub;

      lb = SCIPvarGetLbLocal(mindcons->distvar);
      ub = SCIPvarGetUbLocal(mindcons->distvar);

      distbound = getMinimumScalarVariableProduct(lb, ub, mindcons->scalar, mindcons->issquared);
   }
   distbound += mindcons->distval;

   /* the distance will always be non-negative, i.e., the distbound is always statisfied */
   if( SCIPisLE(scip, distbound, 0.0) )
      return SCIP_OKAY;

   /* compute total sum of squared distances */
   for( i = 0; i < nvars; ++i )
      totalsquareddist += pow(maxDistVarPair(vars1[i], vars2[i]), 2);

   /* propagate bounds of variables */
   for( i = 0; i < nvars && !*infeasible; ++i )
   {
      dist = pow(maxDistVarPair(vars1[i], vars2[i]), 2);
      SCIP_CALL( propagateMindConsSimple1D(scip, vars1[i], vars2[i], distbound - (totalsquareddist - dist), infeasible, nred) );
   }

   return SCIP_OKAY;
}

/** propagation method for mind conss */
static
SCIP_DECL_DISTCONSPROP(propMind)
{
   int ntmpred;
   int c;

   assert(propdata != NULL);
   assert(nred != NULL);
   assert(didrun != NULL);

   *nred = 0;
   *didrun = FALSE;

   if( propdata->nmindconss <= 0 )
      return SCIP_OKAY;

   for( c = 0; c < propdata->nmindconss && !(*infeasible); ++c )
   {
      if( !propdata->dopropmind[c] )
         continue;

      SCIP_CALL( propagateMindConsSimple(scip, propdata->mindconss[c], infeasible, &ntmpred) );

      propdata->dopropmind[c] = FALSE;
      *didrun = TRUE;
      *nred += ntmpred;
   }

   return SCIP_OKAY;
}

/** propagates a 1-dimensional component of minfd cons using a simple propagation method */
static
SCIP_RETCODE propagateMinfdConsSimple1D(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable in pair */
   SCIP_Real             coord,              /**< coordinate in pair */
   SCIP_Real             sqdistbound,        /**< squared distance bound for propagation */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility has been detected */
   int*                  nred                /**< pointer to increment by number of found reductions */
   )
{
   SCIP_Real distbound;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(var != NULL);
   assert(infeasible != NULL);
   assert(nred != NULL);

   if( sqdistbound < 0.0 )
      return SCIP_OKAY;

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   distbound = sqrt(sqdistbound);

   if( SCIPisLE(scip, ub, coord) )
   {
      /* var is to the left of coord */
      SCIP_CALL( improveUb(scip, var, coord - distbound, infeasible, nred) );
   }
   else if( SCIPisLE(scip, coord, lb) )
   {
      /* var is to the right of coord */
      SCIP_CALL( improveLb(scip, var, coord + distbound, infeasible, nred) );
   }
   else if( SCIPisGE(scip, ub, coord) && SCIPisLE(scip, lb, coord) )
   {
      /* var's domain contains coord */
      if( SCIPisLT(scip, coord - lb, distbound) )
      {
         SCIP_CALL( improveLb(scip, var, coord + distbound, infeasible, nred) );
      }
      if( !*infeasible && SCIPisLT(scip, ub - coord, distbound) )
      {
         SCIP_CALL( improveUb(scip, var, coord - distbound, infeasible, nred) );
      }
   }

   return SCIP_OKAY;
}

/** propagates a minfd cons using a simple propagation method */
static
SCIP_RETCODE propagateMinfdConsSimple(
   SCIP*                 scip,               /**< SCIP data structure */
   DISTCONS_MINFIXED*    minfdcons,          /**< minfd cons */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility is detected */
   int*                  nred                /**< pointer to store number of found reductions */
   )
{
   SCIP_Real totalsquareddist = 0.0;
   SCIP_Real sqdistbound;
   SCIP_Real tmpmaxdist;
   SCIP_VAR** vars;
   SCIP_Real* point;
   SCIP_Real dist;
   int nvars;
   int i;

   assert(minfdcons != NULL);
   assert(infeasible != NULL);
   assert(nred != NULL);

   *nred = 0;
   *infeasible = FALSE;

   /* get data */
   vars = minfdcons->vars;
   point = minfdcons->point;
   nvars = minfdcons->lenvars;
   sqdistbound = minfdcons->distval;
   if( minfdcons->hasvardist )
      sqdistbound += getMinimumScalarVariableProduct(SCIPvarGetLbLocal(minfdcons->distvar),
         SCIPvarGetUbLocal(minfdcons->distvar), minfdcons->scalar, FALSE);

   /* do nothing if the bound is negative (this can happen for variable distances if no solution has been found yet) */
   if( SCIPisLE(scip, sqdistbound, 0.0) )
      return SCIP_OKAY;

   /* compute total sum of squared distances */
   for( i = 0; i < nvars; ++i )
   {
      tmpmaxdist = MAX(ABS(SCIPvarGetLbLocal(vars[i]) - point[i]), ABS(SCIPvarGetUbLocal(vars[i]) - point[i]));
      totalsquareddist += pow(tmpmaxdist, 2);
   }

   /* propagate bounds of variables */
   for( i = 0; i < nvars && !*infeasible; ++i )
   {
      dist = MAX(ABS(SCIPvarGetLbLocal(vars[i]) - point[i]), ABS(SCIPvarGetUbLocal(vars[i]) - point[i]));
      dist = SQR(dist);

      SCIP_CALL( propagateMinfdConsSimple1D(scip, vars[i], point[i], sqdistbound - (totalsquareddist - dist),
            infeasible, nred) );
   }

   return SCIP_OKAY;
}

/** propagation method for minfd conss */
static
SCIP_DECL_DISTCONSPROP(propMinfd)
{
   int ntmpred;
   int c;

   assert(propdata != NULL);
   assert(nred != NULL);
   assert(didrun != NULL);

   *nred = 0;
   *didrun = FALSE;

   if( propdata->nminfdconss <= 0 )
      return SCIP_OKAY;

   for( c = 0; c < propdata->nminfdconss && !(*infeasible); ++c )
   {
      if( !propdata->dopropminfd[c] )
         continue;

      SCIP_CALL( propagateMinfdConsSimple(scip, propdata->minfdconss[c], infeasible, &ntmpred) );

      propdata->dopropminfd[c] = FALSE;
      *didrun = TRUE;
      *nred += ntmpred;
   }

   return SCIP_OKAY;
}

/** returns if the intersections on the edges of a box and a ball centered at a point is close enough to a point set */
static
SCIP_Bool intersectionPointsEdgeCloseEnough(
   SCIP*                 scip,               /**< SCIP data structure */
   BOX*                  box,                /**< box for whose edges the intersections are computed */
   SCIP_Real*            point,              /**< point based on which we compute the intersection points */
   SCIP_Real             sqrad,              /**< squared radius of ball centered at point */
   BOX*                  pointset,           /**< point set (vertices of a box) */
   SCIP_Real             sqdist,             /**< squared distance to point set */
   int                   freeidx,            /**< index of dimension which varies along the edge */
   int                   edgeidx,            /**< identifier of lbs/ubs for fixed values along edge */
   int                   dim                 /**< dimension of point */
   )
{
   SCIP_Real intersection[3];
   SCIP_Real* vertex;
   SCIP_Real sqdistval;
   SCIP_Real distval;
   SCIP_Real diff;
   SCIP_Real edgecoord1;
   SCIP_Real edgecoord2;
   int tmpidx;
   int vidx;
   int d;

   assert(scip != NULL);
   assert(box != NULL);
   assert(point != NULL);
   assert(sqrad >= 0);
   assert(pointset != NULL);
   assert(sqdist >= 0);
   assert(0 <= freeidx && freeidx < dim );
   assert(0 <= edgeidx);
   assert(dim != 2 || edgeidx < 2);
   assert(dim != 3 || edgeidx < 4);

   /*
    * compute intersection of affine hull of edge with the ball centered at point
    */
   sqdistval = sqrad;

   /* determine intersection for all coordinates except for the free coordinate of the edge */
   tmpidx = edgeidx;
   for( d = 0; d < dim; ++d )
   {
      SCIP_Real edgecoord;

      /* skip the freeidx, this will be computed later */
      if( d == freeidx )
         continue;

      edgecoord = (tmpidx % 2) == 0 ? box->lbs[d] : box->ubs[d];
      intersection[d] = edgecoord;
      diff = point[d] - edgecoord;
      sqdistval -= diff * diff;
      tmpidx = tmpidx >> 1;
   }

   /*
    * freeidx coordinate is solution of (x - point[freeidx])^2 = sqdistval
    */

   /* there is no intersection or the ball is tangential to the edge
    *
    * The checks outside this function make sure that the other balls contain the endpoints of the edge.
    * This is checked by the DEBUGCHECK.
    */
   if( SCIPisLE(scip, sqdistval, 0.0) )
      goto DEBUGCHECK;
   assert(sqdistval > 0);

   /* compute the two intersections and check whether the ball contains the edge */
   distval = sqrt(sqdistval);
   edgecoord1 = point[freeidx] + distval;
   edgecoord2 = point[freeidx] - distval;

   if( SCIPisLE(scip, edgecoord2, box->lbs[freeidx]) && SCIPisGE(scip, edgecoord1, box->ubs[freeidx]) )
   {
      /* the ball contains the entire edge */
      return TRUE;
   }
   else if( SCIPisLE(scip, box->lbs[freeidx], edgecoord1) && SCIPisLE(scip, edgecoord1, box->ubs[freeidx]) )
   {
      if( SCIPisLE(scip, box->lbs[freeidx], edgecoord2) && SCIPisLE(scip, edgecoord2, box->ubs[freeidx]) )
      {
         /* the edge contains both intersection points, i.e., the balls of the point set need to contain the edge */
         goto DEBUGCHECK;
      }
      else
         intersection[freeidx] = edgecoord1;
   }
   else if( SCIPisLE(scip, box->lbs[freeidx], edgecoord2) && SCIPisLE(scip, edgecoord2, box->ubs[freeidx]) )
      intersection[freeidx] = edgecoord2;
   else
   {
      /* the edge does not intersect with the ball, i.e., the balls of the point set need to contain the edge */
      goto DEBUGCHECK;
   }

   /* check whether distance of the intersection with the point set is small enough */
   for( vidx = 0; vidx < pointset->nvertices; ++vidx )
   {
      vertex = &pointset->vertices[vidx * dim];
      if( SCIPisGT(scip, getDistL2Sq(vertex, intersection, dim), sqdist) )
         return FALSE;
   }

   return TRUE;

 DEBUGCHECK:
   for( vidx = 0; vidx < pointset->nvertices; ++vidx )
   {
      vertex = &pointset->vertices[vidx * dim];
      intersection[freeidx] = box->lbs[freeidx];
      assert(SCIPisLE(scip, getDistL2Sq(intersection, vertex, dim), sqdist + EPS)
         || SCIPisLE(scip, getDistL2Sq(intersection, point, dim), sqdist + EPS));

      intersection[freeidx] = box->ubs[freeidx];
      assert(SCIPisLE(scip, getDistL2Sq(intersection, vertex, dim), sqdist + EPS)
         || SCIPisLE(scip, getDistL2Sq(intersection, point, dim), sqdist + EPS));
   }

   return TRUE;
}

/** returns whether a box is contained in the pairwise union of balls centered at vertices of two other boxes */
static
SCIP_Bool isBoxContainedUnionVertexBalls(
   SCIP*                 scip,               /**< SCIP data structure */
   BOX*                  commonbox,          /**< container holding vertices of common vars */
   BOX*                  box1,               /**< container holding vertices of first box */
   BOX*                  box2,               /**< container holding vertices of second box */
   SCIP_Real             sqrad1,             /**< squared radius of balls centered at vertices of box1 */
   SCIP_Real             sqrad2              /**< squared radius of balls centered at vertices of box2 */
   )
{
   SCIP_Real* vertex1;
   SCIP_Real* vertex2;
   SCIP_Real* vertex;
   int nedgechecks;
   int posv1;
   int posv2;
   int posv;
   int dim;
   int v1;
   int v2;
   int v;
   int f;

   assert(scip != NULL);
   assert(commonbox != NULL);
   assert(box1 != NULL);
   assert(box2 != NULL);
   assert(commonbox->dim == box1->dim && box1->dim == box2->dim);
   assert(sqrad1 >= 0);
   assert(sqrad2 >= 0);

   dim = box1->dim;
   nedgechecks = (dim == 2) ? 2 : 4;

   /* iterate through vertices of box1 */
   for( v1 = 0; v1 < box1->nvertices; ++v1 )
   {
      posv1 = v1 * dim;
      vertex1 = &box1->vertices[posv1];

      /* find the vertices of commonbox that are too far away from the selected vertex of box1 */
      for( v = 0; v < commonbox->nvertices; ++v )
      {
         posv = v * dim;
         vertex = &commonbox->vertices[posv];

         /* if the vertex is too far away from vertex1, it must be close enough to all vertices of the second box */
         if( SCIPisGT(scip, getDistL2Sq(vertex1, vertex, dim), sqrad1) )
         {
            /* iterate through remaining vertices and check their distance */
            for( v2 = 0; v2 < box2->nvertices; ++v2 )
            {
               posv2 = v2 * dim;
               vertex2 = &box2->vertices[posv2];

               /* if the vertex is also too far away from vertex2, bound is not valid */
               if( SCIPisGT(scip, getDistL2Sq(vertex2, vertex, dim), sqrad2) )
                  return FALSE;
            }
         }
      }
   }

   /* also the intersections of the ball at vertex1 and the edges of commonbox must be contained in
    * the balls centered at the vertices of box2
    */
   for( v1 = 0; v1 < box1->nvertices; ++v1 )
   {
      posv1 = v1 * dim;
      vertex1 = &box1->vertices[posv1];

      /* iterate through free dimension of edge */
      for( f = 0; f < dim; ++f )
      {
         /* check all possible fixed coordinates of edge */
         for( v = 0; v < nedgechecks; ++v )
         {
            if( !intersectionPointsEdgeCloseEnough(scip, commonbox, vertex1, sqrad1, box2, sqrad2, f, v, dim) )
               return FALSE;
         }
      }
   }

   return TRUE;
}

/** propagates a facet of the common variables of a minpd or minfpd cons, i.e., tries to find a better bound
 *  for an inequality of type \f$\pm x_i \leq b\f$ when every point in the part that is cut off from the box
 *  is too close to the vertices of box1 and box2 */
static
SCIP_RETCODE propagatePairedDistConsFacet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            commonvars,         /**< common variables of a paired distance constraint (minpd, minfpd) */
   int                   dim,                /**< dimension of distance constraint's ambient space */
   int                   propidx,            /**< coordinate that is propagated */
   SCIP_Bool             doproplb,           /**< whether a lower bound is propagated */
   BOX*                  commonbox,          /**< container to store vertices of box domain of common variables */
   BOX*                  box1,               /**< container to store vertices of box domain of vars1 of distance cons */
   BOX*                  box2,               /**< container to store vertices of box domain of vars2 of distance cons */
   SCIP_Real             sqrad1,             /**< squared radius of first distance constraint */
   SCIP_Real             sqrad2,             /**< squared radius of second distance constraint */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility is detected */
   int*                  nred,               /**< pointer to store number of found reductions */
   SCIP_Real             multbisection,      /**< multiplier used in bisection to compute point in interval */
   SCIP_Real             gapbisection,       /**< absolute gap between bounds to terminate bisection */
   int                   maxniterbisection   /**< maximum number of refinement iterations of bisection */
   )
{
   int propcoords[4];
   int npropcoords = 0;
   SCIP_Real proplb = -SCIP_REAL_MAX;
   SCIP_Real propub = SCIP_REAL_MAX;
   SCIP_Real boundguess;
   SCIP_Real step;
   int maxiter;
   int cnt = 0;
   int nverts;
   int tmpv;
   int pos;
   int v;
   int d;

   assert(scip != NULL);
   assert(commonvars != NULL);
   assert(2 <= dim && dim <= 3);
   assert(0 <= propidx && propidx < dim);
   assert(commonbox != NULL);
   assert(box1 != NULL);
   assert(box2 != NULL);
   assert(sqrad1 >= 0);
   assert(sqrad2 >= 0);
   assert(infeasible != NULL);
   assert(nred != NULL);

   *nred = 0;
   *infeasible = FALSE;

   /* only propagate if there is sufficient gains to be expected */
   if( commonbox->ubs[propidx] - commonbox->lbs[propidx] < 1e-2 )
      return SCIP_OKAY;

   /* collect information about the box domain of the common variables */
   nverts = 1u << dim;
   for( v = 0, pos = 0; v < nverts; ++v )
   {
      tmpv = v;
      for( d = 0; d < dim; ++d, ++pos )
      {
         commonbox->vertices[pos] =
            (tmpv % 2) == 0 ? SCIPvarGetLbLocal(commonvars[d]) : SCIPvarGetUbLocal(commonvars[d]);

         /* keep track of entries that can be changed due to propagation 
          * to find an improvement on the lower bound, we need to find an box that is infeasible after reducing the upper bound
          * and vice versa for the upper bound
          */
         if( d == propidx )
         {
            if( (doproplb && tmpv % 2 == 1) || (!doproplb && tmpv % 2 == 0) )
               propcoords[npropcoords++] = pos;
         }

         tmpv = tmpv >> 1;
      }
   }
   assert(dim != 2 || npropcoords == 1 || npropcoords == 2);
   assert(dim != 3 || npropcoords == 1 || npropcoords == 4);
   for( d = 0; d < dim; ++d )
   {
      commonbox->lbs[d] = SCIPvarGetLbLocal(commonvars[d]);
      commonbox->ubs[d] = SCIPvarGetUbLocal(commonvars[d]);
   }

   /* try to find better bound via bisection */
   proplb = SCIPvarGetLbLocal(commonvars[propidx]);
   propub = SCIPvarGetUbLocal(commonvars[propidx]);
   maxiter = maxniterbisection > 0 ? maxniterbisection : INT_MAX;

   while( propub - proplb > gapbisection && cnt < maxiter )
   {
      ++cnt;
      step = multbisection * (propub - proplb);
      boundguess = doproplb ? proplb + step : propub - step;

      /* update the vertices and bounds of commonbox according to the guess
       *
       * We use common box to encode the box that we want to cut off based on bisection.
       * For this reason, if we want to propagate a lower bound, we update the upper bound
       * of the box to see whether all parts with a worse lower bound can be cut off. */
      for( v = 0; v < npropcoords; ++v )
         commonbox->vertices[propcoords[v]] = boundguess;
      if( doproplb )
         commonbox->ubs[propidx] = boundguess;
      else
         commonbox->lbs[propidx] = boundguess;

      if( isBoxContainedUnionVertexBalls(scip, commonbox, box1, box2, sqrad1, sqrad2) )
      {
         /* if the box can be cut off, then we might even cut off a larger box,
          * so update the bound we want to propagate */
         if( doproplb )
            proplb = boundguess;
         else
            propub = boundguess;
      }
      else
      {
         /* if the box cannot be cut off, then we might cut off a smaller box,
          * so update the bound we do not want to propagate */
         if( doproplb )
            propub = boundguess;
         else
            proplb = boundguess;
      }
   }

   /* possibly improve a variable bound */
   if( doproplb && proplb > SCIPvarGetLbLocal(commonvars[propidx]) )
   {
      SCIP_CALL( improveLb(scip, commonvars[propidx], proplb, infeasible, nred) );
   }
   else if( !doproplb && propub < SCIPvarGetUbLocal(commonvars[propidx]) )
   {
      SCIP_CALL( improveUb(scip, commonvars[propidx], propub, infeasible, nred) );
   }

   return SCIP_OKAY;
}

/** propagates a minpd cons */
static
SCIP_RETCODE propagateMinpdCons(
   SCIP*                 scip,               /**< SCIP data structure */
   DISTCONS_MINPAIR*     minpdcons,          /**< minpd cons */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility is detected */
   int*                  nred,               /**< pointer to store number of found reductions */
   SCIP_Real             multbisection,      /**< multiplier used in bisection to compute point in interval */
   SCIP_Real             gapbisection,       /**< absolute gap between bounds to terminate bisection */
   int                   maxniterbisection   /**< maximum number of refinement iterations of bisection */
   )
{
   SCIP_Real sqrad1;
   SCIP_Real sqrad2;
   BOX commonbox;
   BOX box1;
   BOX box2;
   int ntmpred;
   int nverts;
   int tmpv;
   int pos;
   int i;
   int d;

   assert(scip != NULL);
   assert(minpdcons != NULL);
   assert(infeasible != NULL);
   assert(nred != NULL);
   assert(2 <= minpdcons->dimension && minpdcons->dimension <= 3);

   *nred = 0;
   *infeasible = FALSE;

   /* prepare box domains for common and distinct variables of pair of mind conss */
   nverts = 1u << minpdcons->dimension;
   commonbox.nvertices = nverts;
   commonbox.dim = minpdcons->dimension;
   box1.nvertices = nverts;
   box1.dim = minpdcons->dimension;
   box2.nvertices = nverts;
   box2.dim = minpdcons->dimension;

   /* get vertices of box domains of distinct variables (they will not change in propagation)
    * (do not set upper/lower bounds, not needed for box1 and box2)
    */
   for( i = 0, pos = 0; i < nverts; ++i )
   {
      tmpv = i;
      for( d = 0; d < minpdcons->dimension; ++d, ++pos )
      {
         box1.vertices[pos] =
            (tmpv % 2) == 0 ? SCIPvarGetLbLocal(minpdcons->vars1[d]) : SCIPvarGetUbLocal(minpdcons->vars1[d]);
         box2.vertices[pos] =
            (tmpv % 2) == 0 ? SCIPvarGetLbLocal(minpdcons->vars2[d]) : SCIPvarGetUbLocal(minpdcons->vars2[d]);
         tmpv = tmpv >> 1;
      }
   }

   /* get squared radii used for propagation */
   sqrad1 = minpdcons->constants[0];
   if( minpdcons->distvars[0] != NULL )
   {
      sqrad1 += getMinimumScalarVariableProduct(SCIPvarGetLbLocal(minpdcons->distvars[0]),
         SCIPvarGetUbLocal(minpdcons->distvars[0]), minpdcons->scalars[0], minpdcons->issquared[0]);
   }

   /* do nothing if the bound is negative (this can happen for variable distances if no solution has been found yet) */
   if( SCIPisLE(scip, sqrad1, 0.0) )
      return SCIP_OKAY;
   assert(sqrad1 >= 0);

   sqrad2 = minpdcons->constants[1];
   if( minpdcons->distvars[1] != NULL )
   {
      sqrad2 += getMinimumScalarVariableProduct(SCIPvarGetLbLocal(minpdcons->distvars[1]),
         SCIPvarGetUbLocal(minpdcons->distvars[1]), minpdcons->scalars[1], minpdcons->issquared[1]);
   }

   /* do nothing if the bound is negative (this can happen for variable distances if no solution has been found yet) */
   if( SCIPisLE(scip, sqrad2, 0.0) )
      return SCIP_OKAY;
   assert(sqrad2 >= 0);

   /* propagate all facets of the box domain of the common variables */
   for( d = 0; d < minpdcons->dimension; ++d )
   {
      for( i = 0; i < 2; ++i )
      {
         SCIP_CALL( propagatePairedDistConsFacet(scip, minpdcons->commonvars, minpdcons->dimension,
               d, i == 0, &commonbox, &box1, &box2, sqrad1, sqrad2,
               infeasible, &ntmpred, multbisection, gapbisection, maxniterbisection) );

         if( *infeasible )
            return SCIP_OKAY;
         *nred += ntmpred;
      }
   }

   return SCIP_OKAY;
}

/** propagation method for minpd conss */
static
SCIP_DECL_DISTCONSPROP(propMinpd)
{
   DISTCONS_MINPAIR* minpdcons;
   SCIP_Longint nodeid;
   int ntmpred;
   int c;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(nred != NULL);
   assert(didrun != NULL);

   *nred = 0;
   *didrun = FALSE;

   if( propdata->nminpdconss <= 0 )
      return SCIP_OKAY;
   nodeid = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   for( c = 0; c < propdata->nminpdconss && !(*infeasible); ++c )
   {
      if( !propdata->dopropminpd[c] )
         continue;

      /* do not propagate the constraint if one constraint of the pair has found a reduction */
      minpdcons = propdata->minpdconss[c];
      if( nodeid == minpdcons->dconss[0]->lastpropnode || nodeid == minpdcons->dconss[1]->lastpropnode )
      {
         propdata->dopropminpd[c] = FALSE;
         continue;
      }
      minpdcons->dconss[0]->lastpropnode = nodeid;
      minpdcons->dconss[1]->lastpropnode = nodeid;

      SCIP_CALL( propagateMinpdCons(scip, propdata->minpdconss[c], infeasible, &ntmpred,
            propdata->multbisection, propdata->gapbisection, propdata->maxniterbisection) );

      propdata->dopropminpd[c] = FALSE;
      *didrun = TRUE;
      *nred += ntmpred;
   }

   return SCIP_OKAY;
}

/** propagates a minfpd cons */
static
SCIP_RETCODE propagateMinfpdCons(
   SCIP*                 scip,               /**< SCIP data structure */
   DISTCONS_MINFPAIR*    minfpdcons,         /**< minfpd cons */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility is detected */
   int*                  nred,               /**< pointer to store number of found reductions */
   SCIP_Real             multbisection,      /**< multiplier used in bisection to compute point in interval */
   SCIP_Real             gapbisection,       /**< absolute gap between bounds to terminate bisection */
   int                   maxniterbisection   /**< maximum number of refinement iterations of bisection */
   )
{
   SCIP_Real sqrad1;
   SCIP_Real sqrad2;
   BOX commonbox;
   BOX box1;
   BOX box2;
   int ntmpred;
   int nverts;
   int i;
   int d;

   assert(scip != NULL);
   assert(minfpdcons != NULL);
   assert(infeasible != NULL);
   assert(nred != NULL);
   assert(2 <= minfpdcons->dimension && minfpdcons->dimension <= 3);

   *nred = 0;
   *infeasible = FALSE;

   /* prepare box domains for common and distinct variables of pair of mind conss */
   nverts = pow(2, minfpdcons->dimension);
   commonbox.nvertices = nverts;
   commonbox.dim = minfpdcons->dimension;
   box1.nvertices = 1;
   box1.dim = minfpdcons->dimension;
   box2.nvertices = 1;
   box2.dim = minfpdcons->dimension;

   /* get fixed points of distance constraints (they will not change in propagation)
    * (do not set upper/lower bounds, not needed for box1 and box2)
    */
   for( d = 0; d < minfpdcons->dimension; ++d)
   {
      box1.vertices[d] = minfpdcons->point1[d];
      box2.vertices[d] = minfpdcons->point2[d];
   }

   /* get squared radii used for propagation */
   sqrad1 = minfpdcons->constants[0];
   if( minfpdcons->distvars[0] != NULL )
   {
      sqrad1 += getMinimumScalarVariableProduct(SCIPvarGetLbLocal(minfpdcons->distvars[0]),
         SCIPvarGetUbLocal(minfpdcons->distvars[0]), minfpdcons->scalars[0], FALSE);
   }

   /* this can happen if the distance variable is unbounded (e.g., no solution found yet) */
   if( SCIPisLE(scip, sqrad1, 0.0) )
      return SCIP_OKAY;

   sqrad2 = minfpdcons->constants[1];
   if( minfpdcons->distvars[1] != NULL )
   {
      sqrad2 += getMinimumScalarVariableProduct(SCIPvarGetLbLocal(minfpdcons->distvars[1]),
         SCIPvarGetUbLocal(minfpdcons->distvars[1]), minfpdcons->scalars[1], FALSE);
   }

   /* this can happen if the distance variable is unbounded (e.g., no solution found yet) */
   if( SCIPisLE(scip, sqrad2, 0.0) )
      return SCIP_OKAY;

   /* propagate all facets of the box domain of the common variables */
   for( d = 0; d < minfpdcons->dimension; ++d )
   {
      for( i = 0; i < 2; ++i )
      {
         SCIP_CALL( propagatePairedDistConsFacet(scip, minfpdcons->vars, minfpdcons->dimension,
               d, i == 0, &commonbox, &box1, &box2, sqrad1, sqrad2,
               infeasible, &ntmpred, multbisection, gapbisection, maxniterbisection) );

         if( *infeasible )
            return SCIP_OKAY;
         *nred += ntmpred;
      }
   }

   return SCIP_OKAY;
}

/** propagation method for minfpd conss */
static
SCIP_DECL_DISTCONSPROP(propMinfpd)
{
   DISTCONS_MINFPAIR* minfpdcons;
   SCIP_Longint nodeid;
   int ntmpred;
   int c;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(nred != NULL);
   assert(didrun != NULL);

   *nred = 0;
   *didrun = FALSE;

   if( propdata->nminfpdconss <= 0 )
      return SCIP_OKAY;
   nodeid = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   for( c = 0; c < propdata->nminfpdconss && !(*infeasible); ++c )
   {
      if( !propdata->dopropminfpd[c] )
         continue;

      /* do not propagate the constraint if one constraint of the pair has found a reduction */
      minfpdcons = propdata->minfpdconss[c];
      if( nodeid == minfpdcons->dconss[0]->lastpropnode || nodeid == minfpdcons->dconss[1]->lastpropnode )
      {
         propdata->dopropminpd[c] = FALSE;
         continue;
      }
      minfpdcons->dconss[0]->lastpropnode = nodeid;
      minfpdcons->dconss[1]->lastpropnode = nodeid;

      SCIP_CALL( propagateMinfpdCons(scip, propdata->minfpdconss[c], infeasible, &ntmpred,
            propdata->multbisection, propdata->gapbisection, propdata->maxniterbisection) );

      propdata->dopropminfpd[c] = FALSE;
      *didrun = TRUE;
      *nred += ntmpred;
   }

   return SCIP_OKAY;
}

/*
 * methods for creating propdata
 */

/** clears propagator data */
static
SCIP_RETCODE clearPropData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   int i;

   assert(propdata != NULL);

   SCIP_CALL( dropVarEvents(scip, propdata->eventhdlr,
         propdata->mindconss, propdata->nmindconss, propdata->eventdatamind,
         propdata->minfdconss, propdata->nminfdconss, propdata->eventdataminfd,
         propdata->minpdconss,  propdata->nminpdconss, propdata->eventdataminpd,
         propdata->minfpdconss,  propdata->nminfpdconss, propdata->eventdataminfpd) );
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->dopropmind, propdata->nmindconss);
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->dopropminfd, propdata->nminfdconss);
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->dopropminpd, propdata->nminpdconss);
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->dopropminfpd, propdata->nminfpdconss);

   SCIPfreeBlockMemoryArrayNull(scip, &propdata->eventdatamind, propdata->nmindconss);
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->eventdataminfd, propdata->nminfdconss);
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->eventdataminpd, propdata->nminpdconss);
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->eventdataminfpd, propdata->nminfpdconss);

   for( i = 0; i < propdata->nmindconss; ++i )
   {
      SCIP_CALL( freeDistconsMind(scip, &propdata->mindconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->mindconss, propdata->lenmindconss);
   propdata->nmindconss = 0;
   propdata->lenmindconss = 0;

   for( i = 0; i < propdata->nminfdconss; ++i )
   {
      SCIP_CALL( freeDistconsMinfd(scip, &propdata->minfdconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->minfdconss, propdata->lenminfdconss);
   propdata->nminfdconss = 0;
   propdata->lenminfdconss = 0;

   for( i = 0; i < propdata->nminpdconss; ++i )
   {
      SCIPfreeBlockMemory(scip, &propdata->minpdconss[i]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->minpdconss, propdata->lenminpdconss);
   propdata->nminpdconss = 0;
   propdata->lenminpdconss = 0;

   for( i = 0; i < propdata->nminfpdconss; ++i )
   {
      SCIPfreeBlockMemory(scip, &propdata->minfpdconss[i]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->minfpdconss, propdata->lenminfpdconss);
   propdata->nminfpdconss = 0;
   propdata->lenminfpdconss = 0;

   for( i = 0; i < propdata->nprops; ++i )
   {
      SCIPfreeBlockMemory(scip, &propdata->props[i]);
   }
   propdata->nprops = 0;

   return SCIP_OKAY;
}

/** initializes propagator data */
static
SCIP_RETCODE initPropData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(propdata != NULL);

   BMSclearMemory(propdata);

   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &propdata->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecDistcons, (SCIP_EVENTHDLRDATA*)propdata) );

   return SCIP_OKAY;
}

/*
 * methods for detecting new types of distance constraints from existing distance constraints
 */

/** given the data of two 2-dimensional mind conss, tries to extract common and distinct variables */
static
void tryExtractDataMinpd2D(
   SCIP_VAR**            vars1a,             /**< first variable array of first mind cons */
   SCIP_VAR**            vars2a,             /**< second variable array of first mind cons */
   SCIP_VAR**            vars1b,             /**< first variable array of second mind cons */
   SCIP_VAR**            vars2b,             /**< second variable array of second mind cons */
   SCIP_VAR**            commonvars,         /**< array to store common variables */
   SCIP_VAR**            distinctvars1,      /**< array to store distinct variables of first cons */
   SCIP_VAR**            distinctvars2,      /**< array to store distinct variables of second cons */
   SCIP_Bool*            success             /**< pointer to store whether variables could be extracted successfully */
   )
{
   int i;
   int j;

   *success = TRUE;

   for( i = 0; i < 2; ++i )
   {
      for( j = 0; j < 2; ++j )
      {
         /* check all possible pairs of common variables */
         if( vars1a[i] == vars1b[j] && vars1a[1 - i] == vars1b[1 - j] )
         {
            commonvars[0] = vars1a[i];
            commonvars[1] = vars1a[1 - i];
            distinctvars1[0] = vars2a[i];
            distinctvars1[1] = vars2a[1 - i];
            distinctvars2[0] = vars2b[j];
            distinctvars2[1] = vars2b[1 - j];
         }
         else if( vars1a[i] == vars2b[j] && vars1a[1 - i] == vars2b[1 - j] )
         {
            commonvars[0] = vars1a[i];
            commonvars[1] = vars1a[1 - i];
            distinctvars1[0] = vars2a[i];
            distinctvars1[1] = vars2a[1 - i];
            distinctvars2[0] = vars1b[j];
            distinctvars2[1] = vars1b[1 - j];
         }
         else if( vars2a[i] == vars1b[j] && vars2a[1 - i] == vars1b[1 - j] )
         {
            commonvars[0] = vars2a[i];
            commonvars[1] = vars2a[1 - i];
            distinctvars1[0] = vars1a[i];
            distinctvars1[1] = vars1a[1 - i];
            distinctvars2[0] = vars2b[j];
            distinctvars2[1] = vars2b[1 - j];
         }
         else if( vars2a[i] == vars2b[j] && vars2a[1 - i] == vars2b[1 - j] )
         {
            commonvars[0] = vars2a[i];
            commonvars[1] = vars2a[1 - i];
            distinctvars1[0] = vars1a[i];
            distinctvars1[1] = vars1a[1 - i];
            distinctvars2[0] = vars1b[j];
            distinctvars2[1] = vars1b[1 - j];
         }
      }
   }
   *success = FALSE;
}


/** given the variables of two 2-dimensional minfd cons, tries to extract common and distinct variables */
static
void tryExtractDataMinfpd2D(
   SCIP_VAR**            dcvars1,            /**< variable array of first minfd cons */
   SCIP_Real*            dcpoint1,           /**< point of first minfd cons */
   SCIP_VAR**            dcvars2,            /**< variable array of second minfd cons */
   SCIP_Real*            dcpoint2,           /**< point of second minfd constraint */
   SCIP_VAR**            vars,               /**< array to store common variables */
   SCIP_Real*            point1,             /**< array to store point of first minfd cons */
   SCIP_Real*            point2,             /**< array to store point of second minfd cons */
   SCIP_Bool*            success             /**< pointer to store whether variables could be extracted successfully */
   )
{
   *success = TRUE;

   if( dcvars1[0] == dcvars2[0] && dcvars1[1] == dcvars2[1] )
   {
      vars[0] = dcvars1[0];
      vars[1] = dcvars1[1];
      point1[0] = dcpoint1[0];
      point1[1] = dcpoint1[1];
      point2[0] = dcpoint2[0];
      point2[1] = dcpoint2[1];
   }
   else if( dcvars1[0] == dcvars2[1] && dcvars1[1] == dcvars2[0] )
   {
      vars[0] = dcvars1[0];
      vars[1] = dcvars1[1];
      point1[0] = dcpoint1[0];
      point1[1] = dcpoint1[1];
      point2[0] = dcpoint2[1];
      point2[1] = dcpoint2[0];
   }
   *success = FALSE;
}


/** given the variables of two 3-dimensional mind conss, tries to extract common and distinct variables */
static
void tryExtractDataMinpd3D(
   SCIP_VAR**            vars1a,             /**< first variable array of first mind cons */
   SCIP_VAR**            vars2a,             /**< second variable array of first mind cons */
   SCIP_VAR**            vars1b,             /**< first variable array of second mind cons */
   SCIP_VAR**            vars2b,             /**< second variable array of second mind cons */
   SCIP_VAR**            commonvars,         /**< array to store common variables */
   SCIP_VAR**            distinctvars1,      /**< array to store distinct variables of first mind cons */
   SCIP_VAR**            distinctvars2,      /**< array to store distinct variables of second mind cons */
   SCIP_Bool*            success             /**< pointer to store whether variables could be extracted successfully */
   )
{
   SCIP_VAR** varsdc1[2];
   SCIP_VAR** varsdc2[2];
   int j;
   int offset;
   int ja;
   int jb;
   int r1a;
   int r2a;
   int r3a;
   int r1b;
   int r2b;
   int r3b;

   *success = TRUE;

   /* initialize variable matrices of the two constraints */
   varsdc1[0] = vars1a;
   varsdc1[1] = vars2a;
   varsdc2[0] = vars1b;
   varsdc2[1] = vars2b;

   /* check for all reorderings of variables of second constraint whether they match the variables of the first */
   for( j = 0; j < 3; ++j )
   {
      for( offset = 1; offset < 3; ++offset )
      {
         ja = (j + offset) % 3;
         jb = 3 - j - ja;
         /* here, standard ordering (0,1,2) becomes (j,ja,jb) */

         /* matching variables do not need to be contained in the same subarray vars1/vars2
          * -> check all reorderings of vars1 and vars2 in first and second constraint
          */
         for( r1a = 0; r1a < 2; ++r1a )
         {
            for( r1b = 0; r1b < 2; ++r1b )
            {
               if( varsdc1[r1a][0] != varsdc2[r1b][j] )
                  continue;

               for( r2a = 0; r2a < 2; ++r2a )
               {
                  for( r2b = 0; r2b < 2; ++r2b )
                  {
                     if( varsdc1[r2a][1] != varsdc2[r2b][ja] )
                        continue;

                     for( r3a = 0; r3a < 2; ++r3a )
                     {
                        for( r3b = 0; r3b < 2; ++r3b )
                        {
                           if( varsdc1[r3a][2] != varsdc2[r3b][jb] )
                              continue;

                           /* we have found all common variables */
                           commonvars[0] = varsdc1[r1a][0];
                           commonvars[1] = varsdc1[r2a][1];
                           commonvars[2] = varsdc1[r3a][2];
                           distinctvars1[0] = varsdc1[1 - r1a][0];
                           distinctvars1[1] = varsdc1[1 - r2a][1];
                           distinctvars1[2] = varsdc1[1 - r3a][2];
                           distinctvars2[0] = varsdc2[1 - r1b][0];
                           distinctvars2[1] = varsdc2[1 - r2b][1];
                           distinctvars2[2] = varsdc2[1 - r3b][2];
                        }
                     }
                  }
               }
            }
         }
      }
   }
   *success = FALSE;
}

/** given the data of two 3-dimensional minfd conss, tries to extract common and distinct variables */
static
void tryExtractDataMinfpd3D(
   SCIP_VAR**            dcvars1,            /**< variable array of first minfd cons */
   SCIP_Real*            dcpoint1,           /**< point of first minfd cons */
   SCIP_VAR**            dcvars2,            /**< variable array of second minfd cons */
   SCIP_Real*            dcpoint2,           /**< point of second minfd cons */
   SCIP_VAR**            vars,               /**< array to store common variables */
   SCIP_Real*            point1,             /**< array to store point of first minfd cons */
   SCIP_Real*            point2,             /**< array to store point of second minfd cons */
   SCIP_Bool*            success             /**< pointer to store whether variables could be extracted successfully */
   )
{
   int i;
   int j;
   int k;

   *success = TRUE;

   for( i = 0; i < 3; ++i )
   {
      for( j = 0; j < 3; ++j )
      {
         if( i == j )
            continue;

         k = 3 - i - j;
         if( dcvars1[0] == dcvars2[i] && dcvars1[1] == dcvars2[j] && dcvars1[2] == dcvars2[k] )
         {
            vars[0] = dcvars1[i];
            vars[1] = dcvars1[j];
            vars[2] = dcvars1[k];
            point1[0] = dcpoint1[i];
            point1[1] = dcpoint1[j];
            point1[2] = dcpoint1[k];
            point2[0] = dcpoint2[i];
            point2[1] = dcpoint2[j];
            point2[2] = dcpoint2[k];
         }
      }
   }
   *success = FALSE;
}

/** tries whether pairs of mind conss can be propagated together */
static
SCIP_RETCODE tryCreateMinpdCons(
   SCIP*                 scip,               /**< SCIP data structure */
   DISTCONS_MINPAIR***   minpdconss,         /**< pointer to mind conss */
   int*                  nminpdconss,        /**< pointer to store number of mind conss */
   int*                  lenminpdconss,      /**< pointer to store length of *minpdconss */
   DISTCONS_MIN*         dc1,                /**< first mind cons to be paired */
   DISTCONS_MIN*         dc2                 /**< second mind cons to be paired */
   )
{
   DISTCONS_MINPAIR* mdc;
   SCIP_VAR* commonvars[3];
   SCIP_VAR* vars1[3];
   SCIP_VAR* vars2[3];
   SCIP_Bool success;
   int dim;
   int d;

   assert(scip != NULL);
   assert(minpdconss != NULL);
   assert(nminpdconss != NULL);
   assert(lenminpdconss != NULL);
   assert(dc1 != NULL);
   assert(dc2 != NULL);
   assert(dc1->lenvars == dc2->lenvars);
   assert(dc1->lenvars == 2 || dc1->lenvars == 3);

   dim = dc1->lenvars;

   if( dim == 2 )
   {
      tryExtractDataMinpd2D(dc1->vars1, dc1->vars2, dc2->vars1, dc2->vars2, commonvars, vars1, vars2, &success);
   }
   else
   {
      assert(dim == 3);

      tryExtractDataMinpd3D(dc1->vars1, dc1->vars2, dc2->vars1, dc2->vars2, commonvars, vars1, vars2, &success);
   }

   if( !success )
      return SCIP_OKAY;

   /* possibly reallocate memory */
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, minpdconss, lenminpdconss, *nminpdconss + 1) );

   /* store information about pair of mind conss */
   SCIP_CALL( SCIPallocBlockMemory(scip, &(*minpdconss)[*nminpdconss]) );
   mdc = (*minpdconss)[*nminpdconss];
   mdc->dimension = dim;
   for( d = 0; d < dim; ++d )
      mdc->commonvars[d] = commonvars[d];
   for( d = 0; d < dim; ++d )
      mdc->vars1[d] = vars1[d];
   for( d = 0; d < dim; ++d )
      mdc->vars2[d] = vars2[d];
   mdc->distvars[0] = dc1->distvar;
   mdc->distvars[1] = dc2->distvar;
   mdc->scalars[0] = dc1->scalar;
   mdc->scalars[1] = dc2->scalar;
   mdc->constants[0] = dc1->distval;
   mdc->constants[1] = dc2->distval;
   mdc->issquared[0] = dc1->issquared;
   mdc->issquared[1] = dc2->issquared;
   mdc->dconss[0] = dc1;
   mdc->dconss[1] = dc2;
   ++*nminpdconss;

   return SCIP_OKAY;
}

/** tries whether pairs of minfd conss can be propagated together */
static
SCIP_RETCODE tryCreateMinfpdCons(
   SCIP*                 scip,               /**< SCIP data structure */
   DISTCONS_MINFPAIR***  minfpdconss,        /**< pointer to minfpd conss */
   int*                  nminfpdconss,       /**< pointer to store number of minfpd conss */
   int*                  lenminfpdconss,     /**< pointer to store length of *minfpdconss */
   DISTCONS_MINFIXED*    dc1,                /**< first minfd cons to be paired */
   DISTCONS_MINFIXED*    dc2                 /**< second minfd cons to be paired */
   )
{
   DISTCONS_MINFPAIR* mdc;
   SCIP_VAR* vars[3];
   SCIP_Real point1[3];
   SCIP_Real point2[3];
   SCIP_Bool success;
   int dim;
   int d;

   assert(scip != NULL);
   assert(minfpdconss != NULL);
   assert(nminfpdconss != NULL);
   assert(lenminfpdconss != NULL);
   assert(dc1 != NULL);
   assert(dc2 != NULL);
   assert(dc1->lenvars == dc2->lenvars);
   assert(dc1->lenvars == 2 || dc1->lenvars == 3);

   dim = dc1->lenvars;

   if( dim == 2 )
   {
      tryExtractDataMinfpd2D(dc1->vars, dc1->point, dc2->vars, dc2->point, vars, point1, point2, &success);
   }
   else
   {
      assert(dim == 3);

      tryExtractDataMinfpd3D(dc1->vars, dc1->point, dc2->vars, dc2->point, vars, point1, point2, &success);
   }

   if( !success )
      return SCIP_OKAY;

   /* possibly reallocate memory */
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, minfpdconss, lenminfpdconss, *nminfpdconss + 1) );

   /* store information about minfpd cons */
   SCIP_CALL( SCIPallocBlockMemory(scip, &(*minfpdconss)[*nminfpdconss]) );
   mdc = (*minfpdconss)[*nminfpdconss];
   mdc->dimension = dim;
   for( d = 0; d < dim; ++d )
      mdc->vars[d] = vars[d];
   for( d = 0; d < dim; ++d )
      mdc->point1[d] = point1[d];
   for( d = 0; d < dim; ++d )
      mdc->point2[d] = point2[d];
   mdc->distvars[0] = dc1->distvar;
   mdc->distvars[1] = dc2->distvar;
   mdc->scalars[0] = dc1->scalar;
   mdc->scalars[1] = dc2->scalar;
   mdc->constants[0] = dc1->distval;
   mdc->constants[1] = dc2->distval;
   mdc->dconss[0] = dc1;
   mdc->dconss[1] = dc2;
   ++*nminfpdconss;

   return SCIP_OKAY;
}

/** detects minpd conss */
static
SCIP_RETCODE detectMinpdConss(
   SCIP*                 scip,               /**< SCIP data structure */
   DISTCONS_MIN**        mindconss,          /**< mind conss */
   int                   nmindconss,         /**< number of mind conss */
   DISTCONS_MINPAIR***   minpdconss,         /**< pointer to hold minpd conss */
   int*                  nminpdconss,        /**< pointer to hold number of minpd conss */
   int*                  lenminpdconss       /**< pointer to store length of *minpdconss */
   )
{
   DISTCONS_MIN* dc1;
   DISTCONS_MIN* dc2;
   int c1;
   int c2;
   int dim;

   assert(scip != NULL);
   assert(mindconss != NULL || nmindconss == 0);
   assert(minpdconss != NULL);
   assert(nminpdconss != NULL);
   assert(lenminpdconss != NULL);

   /* if there are not enough mind, we cannot pair them */
   if( nmindconss < 2 )
      return SCIP_OKAY;

   /* iterate through tuples of mind conss and try to pair them */
   for( c1 = 0; c1 < nmindconss; ++c1 )
   {
      dc1 = mindconss[c1];
      assert(dc1 != NULL);

      dim = dc1->lenvars;
      assert(dim >= 1);

      /* currently, we only support dimension 2 or 3 */
      if( dim < 2 || dim > 3 )
         continue;

      for( c2 = c1 + 1; c2 < nmindconss; ++c2 )
      {
         dc2 = mindconss[c2];
         assert(dc2 != NULL);

         /* if the dimension of the constraints' ambient spaces do not match, we cannot pair them */
         if( dim != dc2->lenvars )
            continue;

         SCIP_CALL( tryCreateMinpdCons(scip, minpdconss, nminpdconss, lenminpdconss, dc1, dc2) );
      }
   }

   return SCIP_OKAY;
}

/** detects minfpd conss */
static
SCIP_RETCODE detectMinfdpConss(
   SCIP*                 scip,               /**< SCIP data structure */
   DISTCONS_MINFIXED**   minfdconss,         /**< minfd conss */
   int                   nminfdconss,        /**< number of minfd conss */
   DISTCONS_MINFPAIR***  minfpdconss,        /**< pointer to hold minfpd conss */
   int*                  nminfpdconss,       /**< pointer to hold number of minfpd conss */
   int*                  lenminfpdconss      /**< pointer to store length of *minfdpconss */
   )
{
   DISTCONS_MINFIXED* dc1;
   DISTCONS_MINFIXED* dc2;
   int c1;
   int c2;
   int dim;

   assert(scip != NULL);
   assert(minfdconss != NULL || nminfdconss == 0);
   assert(minfpdconss != NULL);
   assert(nminfpdconss != NULL);
   assert(lenminfpdconss != NULL);

   /* if there are not enough minfd conss, we cannot pair them */
   if( nminfdconss < 2 )
      return SCIP_OKAY;

   /* iterate through tuples of minfd conss and try to pair them */
   for( c1 = 0; c1 < nminfdconss; ++c1 )
   {
      dc1 = minfdconss[c1];
      assert(dc1 != NULL);

      dim = dc1->lenvars;
      assert(dim >= 1);

      /* currently, we only support dimension 2 or 3 */
      if( dim < 2 || dim > 3 )
         continue;

      for( c2 = c1 + 1; c2 < nminfdconss; ++c2 )
      {
         dc2 = minfdconss[c2];
         assert(dc2 != NULL);

         /* if the dimension of the constraints' ambient spaces do not match, we cannot pair them */
         if( dim != dc2->lenvars )
            continue;

         SCIP_CALL( tryCreateMinfpdCons(scip, minfpdconss, nminfpdconss, lenminfpdconss, dc1, dc2) );
      }
   }

   return SCIP_OKAY;
}

/** registers active distance propagators */
static
SCIP_RETCODE registerPropagators(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of distance propagator */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);

   if( propdata->nmindconss > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &propdata->props[propdata->nprops]) );
      propdata->props[propdata->nprops++]->prop = propMind;
   }
   if( propdata->nminfdconss > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &propdata->props[propdata->nprops]) );
      propdata->props[propdata->nprops++]->prop = propMinfd;
   }
   if( propdata->nminpdconss > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &propdata->props[propdata->nprops]) );
      propdata->props[propdata->nprops++]->prop = propMinpd;
   }
   if( propdata->nminfpdconss > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &propdata->props[propdata->nprops]) );
      propdata->props[propdata->nprops++]->prop = propMinfpd;
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of propagator
 */

/** deinitialization method of propagator (called before transformed problem is freed) */
static
SCIP_DECL_PROPEXIT(propExitDistance)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIP_CALL( clearPropData(scip, propdata) );

   return SCIP_OKAY;
}

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeDistance)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIP_CALL( clearPropData(scip, propdata) );

   SCIPfreeBlockMemory(scip, &propdata);

   return SCIP_OKAY;
}

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecDistance)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Bool didrun = FALSE;
   SCIP_Bool infeasible = FALSE;
   SCIP_Bool tmpdidrun;
   int ntmpred = 0;
   int nred = 0;
   int p;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   *result = SCIP_DIDNOTRUN;

   if( !propdata->detectedconss )
   {
      int maxndconss;

      maxndconss = propdata->maxndconssforpairs;

      /* detect distance constraints */
      SCIP_CALL( getMinDistanceConss(scip, &propdata->mindconss, &propdata->nmindconss,
            &propdata->lenmindconss, &propdata->minfdconss, &propdata->nminfdconss, &propdata->lenminfdconss) );

      if( maxndconss > 0 && maxndconss >= propdata->nmindconss + propdata->nminfdconss )
      {
         SCIP_CALL( detectMinpdConss(scip, propdata->mindconss, propdata->nmindconss,
               &propdata->minpdconss, &propdata->nminpdconss, &propdata->lenminpdconss) );

         SCIP_CALL( detectMinfdpConss(scip, propdata->minfdconss, propdata->nminfdconss,
               &propdata->minfpdconss, &propdata->nminfpdconss, &propdata->lenminfpdconss) );
      }
      propdata->detectedconss = TRUE;

      SCIP_CALL( catchVarEvents(scip, propdata->eventhdlr,
            propdata->mindconss, propdata->nmindconss, &propdata->dopropmind, &propdata->eventdatamind,
            propdata->minfdconss, propdata->nminfdconss, &propdata->dopropminfd, &propdata->eventdataminfd,
            propdata->minpdconss,  propdata->nminpdconss, &propdata->dopropminpd, &propdata->eventdataminpd,
            propdata->minfpdconss,  propdata->nminfpdconss, &propdata->dopropminfpd, &propdata->eventdataminfpd) );

      /* register active propagators */
      SCIP_CALL( registerPropagators(scip, propdata) );
   }

   /* call the different propagation methods */
   for( p = 0; p < propdata->nprops; ++p )
   {
      SCIP_CALL( propdata->props[p]->prop(scip, propdata, &tmpdidrun, &infeasible, &ntmpred) );

      if( infeasible )
         break;
      nred += ntmpred;
      didrun = didrun || tmpdidrun;
   }

   if( infeasible )
      *result = SCIP_CUTOFF;
   else if( nred > 0 )
      *result = SCIP_REDUCEDDOM;
   else if( didrun )
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/*
 * propagator specific interface methods
 */

/** creates the distance propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropDistance(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create propagator data */
   propdata = NULL;
   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );

   SCIP_CALL( initPropData(scip, propdata) );

   prop = NULL;

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecDistance, propdata) );

   assert(prop != NULL);

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropExit(scip, prop, propExitDistance) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeDistance) );

   SCIP_CALL( SCIPaddIntParam(scip, "propagating/" PROP_NAME "/maxndconssforpairs",
         "maximum number of distance constraints for which propagation of pairs is applied",
         &propdata->maxndconssforpairs, TRUE, DEFAULT_MAXNDCONSSFORPAIRS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/" PROP_NAME "/multbisection",
         "multiplier used in bisection to split interval [l,u] at point l + multiplier*(u - l)",
         &propdata->multbisection, TRUE, DEFAULT_MULTBISECTION, 0.0001, 0.9999, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/" PROP_NAME "/gapbisection",
         "domain width at which bisection is stopped",
         &propdata->gapbisection, TRUE, DEFAULT_GAPBISECTION, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "propagating/" PROP_NAME "/maxniterbisection",
         "maximum number of iterations used in bisection (0: unbounded)",
         &propdata->maxniterbisection, TRUE, DEFAULT_MAXNITERBISECTION, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
